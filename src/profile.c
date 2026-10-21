/* 
 * This file is part of OpenPMX (https://github.com/deleveld/openpmx).
 * Copyright (c) 2024 Douglas Eleveld.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "openpmx.h"
#include "popmodel.h"
#include "pmxstate.h"
#include "encode.h"
#include "defines.h"
#include "omegafixed.h"
#include "print.h"
#include "githash.h"
#include "utils/c22.h"
#include "utils/errctx.h"
#include "utils/filecontent.h"

#include <gsl/gsl_roots.h>
#include <gsl/gsl_cdf.h>

/// This file implements a function that allows likelihood profiling.
/// It is available in openpmxtran as `profile()`.

/*---------------------------------------------------------------------*/
/* parameter access */
/*---------------------------------------------------------------------*/
static int param_index(const PARAM* const param, const bool _offset1)
{
	return param->index - (_offset1 ? 1 : 0);
}

static const char* param_name(const PARAM* const param)
{
	switch (param->type) {
		case PARAM_THETA:
			return "theta";
		case PARAM_OMEGA:
			return "omega";
		case PARAM_SIGMA:
			return "sigma";
		default:
			return "unknown";
	}
}

static int param_count(const POPMODEL* const popmodel, const PARAM* const param)
{
	switch (param->type) {
		case PARAM_THETA:
			return popmodel->ntheta;
		case PARAM_OMEGA:
			return popmodel->nomega;
		case PARAM_SIGMA:
			return popmodel->nsigma;
		default:
			return 0;
	}
}

static double* param_value(POPMODEL* const popmodel, const PARAM* const param, const int i)
{
	switch (param->type) {
		case PARAM_THETA:
			return &popmodel->theta[i];
		case PARAM_OMEGA:
			return &popmodel->omega[i][i];
		case PARAM_SIGMA:
			return &popmodel->sigma[i];
		default:
			return NULL;
	}
}

/* set the profiled parameter to its test value and fix it */
static void popmodel_apply_param(POPMODEL* const popmodel, const PARAM* const param, const int i)
{
	*param_value(popmodel, param, i) = param->value;
	switch (param->type) {
		case PARAM_THETA:
			popmodel->thetaestim[i] = FIXED;
			break;
		case PARAM_OMEGA:
			popmodel->omegafixed[i][i] = OMEGAFIXED_FIXED;
			break;
		case PARAM_SIGMA:
			popmodel->sigmafixed[i] = 1;
			break;
		default: break;
	}
}

static POPMODEL profile_popmodel_init(const OPENPMX* const source,
									  const PROFILECONFIG* const args,
									  ERRCTX* errctx)
{
	/* first, basic error checking */
	if (args->param.type == PARAM_INVALID) 
		errctx_add(errctx, "%s: profile type invalid\n", __func__);
	if (args->dobjfn < 0.) 
		errctx_add(errctx, "%s: dobjfn cannot be negative\n", __func__);
	if (source->result.type == OBJFN_INVALID ||
		source->result.objfn == DBL_MAX)
		errctx_add(errctx, "%s: model objfn must be valid to do profile\n", __func__);

	/* get popmodel and check for errors */
	var popmodel = popmodel_init(source->theta, source->omega, source->sigma, errctx);
	if (errctx->len)
		return (POPMODEL) { 0 };

	let i = param_index(&args->param, source->data._offset1);
	if (i < 0 || i >= param_count(&popmodel, &args->param)) {
		errctx_add(errctx, "%s: profile %s index (%i) out of bounds\n", 
				   __func__, param_name(&args->param), args->param.index);
		return (POPMODEL) { 0 };
	}
	if (args->param.type == PARAM_OMEGA && popmodel.omegafixed[i][i] == OMEGAFIXED_SAME) {
		errctx_add(errctx, "%s: cannot profile part of omega same block\n", __func__);
		return (POPMODEL) { 0 };
	}
	return popmodel;
}

/*---------------------------------------------------------------------*/
/* iteration file */
/*---------------------------------------------------------------------*/
static FILE* iter_open(const char* sourcename)
{
	if (!sourcename || strlen(sourcename) == 0)
		return 0;
	
	char filename[PATH_MAX];
	snprintf(filename, sizeof(filename), "%s" OPENPMX_PROFILE, sourcename);

	var stream = fopen(filename, "a");
	if (!stream)
		fatal(0, "%s: failed to open \"%s\"\n", __func__, filename);

	/* ftell() on a fresh "a" stream is not portable before the first 
	 * write, so seek explicitly */
	fseek(stream, 0, SEEK_END);
	if (ftell(stream) == 0)
		fprintf(stream, OPENPMX_SFORMAT OPENPMX_SFORMAT 
						OPENPMX_HEADER_FORMAT OPENPMX_HEADER_FORMAT OPENPMX_HEADER_FORMAT 
						OPENPMX_SFORMAT "\n",
						"TYPE", "INDEX", "VALUE", "OBJFN", "TARGET", "NEVAL");
	return stream;
}

static void iter_write(FILE* stream,
					   const PROFILECONFIG* const args,
					   const double objfn, 
					   const double target, 
					   const int neval)
{
	if (!stream)
		return;
	fprintf(stream, OPENPMX_SFORMAT OPENPMX_IFORMAT 
					OPENPMX_TABLE_FORMAT OPENPMX_TABLE_FORMAT OPENPMX_TABLE_FORMAT 
					OPENPMX_IFORMAT "\n",
					param_name(&args->param), args->param.index, args->param.value, 
					objfn, target, neval);
}

/*---------------------------------------------------------------------*/
/* evaluation at a fixed value of the profiled parameter */
/*---------------------------------------------------------------------*/
static void profile_evaluate(OPENPMX* const ret,
							 POPMODEL* const popmodel,
							 const OPENPMX* const source,
							 PROFILECONFIG* const args)
{
	/* never estimate at a non-finite value, callers must check */
	if (!isfinite(args->param.value))
		fatal(0, "%s: profile value %g is not finite\n", __func__, args->param.value);

	/* modify the popmodel */
	popmodel_apply_param(popmodel, &args->param, param_index(&args->param, source->data._offset1));
	memset(&popmodel->result, 0, sizeof(popmodel->result));

	/* reset pmx object, the other parameters start from the source */
	pmx_release_state(ret);
	pmx_copy_popmodel(ret, popmodel);
	pmx_ensure_state(ret);
	if (source->state) {
		/* real checks, asserts vanish with NDEBUG and then we would 
		 * copy etas between mismatched states */
		if (ret->state->idata.nindivid != source->state->idata.nindivid ||
			ret->state->idata.nomega != source->state->idata.nomega)
			fatal(0, "%s: source and profile state sizes differ\n", __func__);

		/* copy the eta state from source to us */
		var etas = idata_etas_alloc(&source->state->idata);
		idata_etas_copy(&etas, &source->state->idata);
		idata_etas_write(&ret->state->idata, &etas);
		idata_etas_free(&etas);
	}

	pmx_estimate(ret, &args->estimate);
}

/*---------------------------------------------------------------------*/
/* root finding */
/*---------------------------------------------------------------------*/
typedef struct {
	const double left;
	const double right;
	const double objfn_target;
	const OPENPMX* const source;
	PROFILECONFIG* args;
	POPMODEL* popmodel;
	OPENPMX* ret;
	int neval;
	bool converged;
	FILE* stream;
	bool final_proposal;
} PROFILEPARAMS;

/* the solver works in normalised x in [0, 1]. sqrt(x) makes the root
 * finding more linear and efficient and, because x is in [0, 1], there
 * is no domain issue. */
static double profile_value(const PROFILEPARAMS* const p, const double x)
{
	return p->left + sqrt(x) * (p->right - p->left);
}

static double root_function(double x, void* _params)
{
	var p = (PROFILEPARAMS*)_params;

	/* get the final proposal without evaluation, so we have a proposed
	 * final estimate which makes more sense as an estimate of the
	 * profile crossing point */
	if (p->final_proposal)
		return 0;

	let args = p->args;
	args->param.value = profile_value(p, x);

	/* x == 0 is the left bracket which is the source fit, no need to
	 * evaluate it again */
	let is_left = (x == 0.);
	var objfn = p->source->result.objfn;
	if (!is_left) {
		profile_evaluate(p->ret, p->popmodel, p->source, args);
		objfn = p->ret->result.objfn;
		p->neval += 1;
	}

	printf("profile eval %s %i value %g objfn %.6f target %.6f neval %i\n", 
		   param_name(&args->param), args->param.index, 
		   args->param.value, objfn, p->objfn_target, p->neval);

	iter_write(p->stream, args, objfn, p->objfn_target, p->neval);

	let err = objfn - p->objfn_target;
	if (!is_left && fabs(err) < args->dobjfn_tol)
		p->converged = true;
	return err;
}

static bool profile_side(const OPENPMX* const source, 
						 PROFILECONFIG* const args, 
						 POPMODEL* const popmodel,
						 OPENPMX* const ret)
{
/// If .dobjfn is 0. then it will be set to the chi-squared distribution 
/// for alpha 0.01 and 1 degree of freedom. This is approximately 6.63.
#define DEFAULT_ALPHA	0.01
#define DEFAULT_DF		1
	if (!(args->dobjfn >= 0.))
		args->dobjfn = gsl_cdf_chisq_Pinv(1.0 - DEFAULT_ALPHA, DEFAULT_DF);
/// The default tolerance is 1 objfn units. This is usually enough that a
/// smoothed line through the profile points is close to the confidence
/// limits.
	if (!(args->dobjfn_tol >= 0.))
		args->dobjfn_tol = 1;
/// The default number of root finding evaluations is 10.
	if (args->maxeval <= 0)
		args->maxeval = 10;

	let left_value = *param_value(popmodel, &args->param, param_index(&args->param, source->data._offset1));
	let right_value = args->param.value;
	if (!isfinite(left_value) || !isfinite(right_value)) {
		/* nothing is estimated and args->param.value is left as given */
		warning(0, "%s: profile value %g or estimate %g not finite, no estimation done\n", 
				__func__, right_value, left_value);
		return false;
	}
	if (left_value == right_value) {
		warning(0, "%s: profile value equals the estimate, nothing to bracket\n", __func__);
		return false;
	}

	var stream = iter_open(source->filename);
	var params = (PROFILEPARAMS) {
		.left = left_value,
		.right = right_value,
		.objfn_target = source->result.objfn + args->dobjfn,
		.source = source,
		.args = args,
		.popmodel = popmodel,
		.ret = ret,
		.stream = stream,
	};
	var F = (gsl_function) {
		.function = root_function,
		.params = &params,
	};

	let type = param_name(&args->param);
	printf("profile start %s %i lower %g upper %g\n"
		   "profile objfn %.6f target %.6f dobjfn %.6f tol %g\n", 
		   type, args->param.index, fmin(left_value, right_value), fmax(left_value, right_value), 
		   source->result.objfn, params.objfn_target, args->dobjfn, args->dobjfn_tol);

	/* gsl_root_fsolver_set() calls root_function at both ends of the
	 * bracket. The left is taken from the source, so only the right is 
	 * really evaluated. */
	var s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
	if (!s)
		fatal(0, "%s: profile root solver alloc failed\n", __func__);

	let oldhandler = gsl_set_error_handler_off();
	let status = gsl_root_fsolver_set(s, &F, 0., 1.);
	if (status != GSL_SUCCESS) {
		let open = (status != GSL_EINVAL) ? NAN
				 : (right_value > left_value) ? INFINITY : -INFINITY;
		printf("profile bound %g failed (%s), limit is %g\n",
			   right_value, gsl_strerror(status), open);
		args->param.value = args->param.lower = args->param.upper = open;
		iter_write(stream, args, 0., params.objfn_target, 0);
		goto cleanup;
	}

	/* iterate the root finder */
	bool solver_ok = true;
	while (!params.converged && params.neval < args->maxeval) {
		let st = gsl_root_fsolver_iterate(s);
		if (st != GSL_SUCCESS && st != GSL_CONTINUE) {
			printf("profile root status no continue (%s)\n", gsl_strerror(st));
			solver_ok = false;
			break;
		}
	}
	if (!params.converged && params.neval >= args->maxeval)
		printf("profile maxeval (%i) reached\n", args->maxeval);

	/* final estimate, but untested into the output file so the users
	 * gets the interpolated best guess so far. This can only be done
	 * by calling gsl_root_fsolver_iterate() which must call the
	 * root_function(). So we add a flag for this special case and we
	 * can avoid doing an expensive evaluation */
	if (solver_ok) {
		params.final_proposal = true;
		gsl_root_fsolver_iterate(s);
	}

	let lower_val	= profile_value(&params, gsl_root_fsolver_x_lower(s));
	let value 		= profile_value(&params, gsl_root_fsolver_root(s));
	let upper_val	= profile_value(&params, gsl_root_fsolver_x_upper(s));
	args->param.lower = fmin(lower_val, upper_val);
	args->param.value = value;
	args->param.upper = fmax(lower_val, upper_val);		

	printf("profile %s value %g [%g, %g]\n", 
			params.converged ? "converged" : "failed",
			args->param.value,
			args->param.lower,
			args->param.upper);
	iter_write(stream, args, 0., params.objfn_target, 0);

cleanup:
	gsl_set_error_handler(oldhandler);
	gsl_root_fsolver_free(s);
	if (stream)
		fclose(stream);

	return true;
}

/// If .value is 0 and .lower < .upper then two profiles are run, one to
/// the .lower bound and one to the .upper bound. On return .lower is the
/// crossing point from the lower profile and .upper is the crossing point
/// from the upper profile. The .value is set to the source estimate. 
/// Otherwise a single profile is run to .value. If a bound does not 
/// cross the .dobjfn threshold then the limit on that side is open and 
/// is returned as -INFINITY (lower) or +INFINITY (upper).
OPENPMX pmx_profile(const OPENPMX* const source, PROFILECONFIG* const args)
{
	ERRCTX errctx = { 0 };
	var popmodel = profile_popmodel_init(source, args, &errctx); 
	if (errctx.len) 
		fatal(0, "%s: %s", __func__, errctx.errmsg);

	var ret = pmx_copy(source);
	if (args->param.value == 0. && args->param.lower < args->param.upper) {
		POPMODEL popmodel_upper = popmodel;
		let value = *param_value(&popmodel, &args->param, param_index(&args->param, source->data._offset1));

		let lo = args->param.lower;
		let hi = args->param.upper;
		
		/* profile the lower bound */
		args->param.lower = 0.;
		args->param.value = lo;
		args->param.upper = 0.;
		profile_side(source, args, &popmodel, &ret);
		let lower_root = args->param.value;

		/* profile the upper bound */
		args->param.lower = 0.;
		args->param.value = hi;
		args->param.upper = 0.;
		profile_side(source, args, &popmodel_upper, &ret);
		let upper_root = args->param.value;

		/* return both bounds */
		args->param.lower = lower_root;
		args->param.value = value;
		args->param.upper = upper_root;
	} else
		profile_side(source, args, &popmodel, &ret);

	ret.filename = 0; /* dont let filename dangle */
	pmx_release_state(&ret);
	return ret;
}

/// pmx_profile_add() evaluates the profile at the single value in args->param.value and
/// append the point to the profile file. The profiled parameter is fixed
/// and the others are re-estimated starting from the source fit. Only
/// args->param and args->estimate are used.
OPENPMX pmx_profile_add(const OPENPMX* const source, PROFILECONFIG* const args)
{
	ERRCTX errctx = { 0 };
	var popmodel = profile_popmodel_init(source, args, &errctx);
	if (errctx.len)
		fatal(0, "%s: %s", __func__, errctx.errmsg);

	if (!isfinite(args->param.value)) {
		warning(0, "%s: profile value %g not finite, no estimation done\n", 
				__func__, args->param.value);
		return (OPENPMX){ 0 };
	}

	var ret = pmx_copy(source);
	FILE* stream = iter_open(source->filename);

	printf("profile add %s %i value %g\n",
		   param_name(&args->param), args->param.index, args->param.value);
	profile_evaluate(&ret, &popmodel, source, args);
	let objfn = ret.result.objfn;
	printf("profile add objfn %.6f\n", objfn);
	iter_write(stream, args, objfn, 0., 1);

	if (stream)
		fclose(stream);
	ret.filename = 0; /* dont let filename dangle */
	pmx_release_state(&ret);
	return ret;
}

void pmx_profile_clean(const OPENPMX* const openpmx, const PARAM* const param)
{
	if (!openpmx || !openpmx->filename || strlen(openpmx->filename) == 0)
		return;

	char filename[PATH_MAX];
	snprintf(filename, sizeof(filename), "%s" OPENPMX_PROFILE, openpmx->filename);

	/* If type passed is PARAM_INVALID, delete the whole .profile.iter file */
	if (param->type == PARAM_INVALID) {
		remove(filename);
		return;
	}

	/* read in the .profile.iter file */
	ERRCTX errctx = { 0 };
	var itertext = filecontent_alloc(filename, &errctx);
	if (errctx.len)
		fatal(0, "%s: %s", __func__, errctx.errmsg);

	/* we dont need a temp file, we can directly re-write the .profile.iter */
	var fout = fopen(filename, "w");
	if (!fout) 
		goto cleanup;

	let target_name = param_name(param);
	let target_index = param->index; /* we dont have to adjust for _offset */
	var ptr = itertext.mutptr;
	while (*ptr != '\0') {
		// Find the end of the current line
		var next_line = strchr(ptr, '\n');
		if (next_line) 
			*next_line = '\0'; // Temporarily terminate the line

		char line_type[64];
		int line_index = -1;
		let scanned = sscanf(ptr, "%63s %d", line_type, &line_index);
		let match = (scanned == 2 && strcmp(line_type, target_name) == 0 && line_index == target_index);

		// Write out non-matching lines
		if (!match && *ptr != '\0') 
			fprintf(fout, "%s\n", ptr);

		if (!next_line) 
			break;
		ptr = next_line + 1;
	}
	fclose(fout);
cleanup:
	filecontent_free(&itertext);
}
