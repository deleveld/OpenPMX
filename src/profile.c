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

#include <limits.h>
#include <string.h>
#include <assert.h>

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

#include <gsl/gsl_roots.h>
#include <gsl/gsl_cdf.h>

/// This file implements a function that allows likelihood profiling.
/// It is available in openpmxtran as `profile()`.

static POPMODEL profile_popmodel_init(const OPENPMX* const source,
									  PROFILECONFIG* const args,
									  ERRCTX* errctx)
{
	/* first, basic error checking */
	if (args->type == PROFILE_INVALID) 
		errctx_add(errctx, "%s: profile type invalid\n", __func__);
	if (args->dobjfn < 0.) 
		errctx_add(errctx, "%s: dobjfn must be positive\n", __func__);
	if (source->result.type == OBJFN_INVALID ||
		source->result.objfn == DBL_MAX)
		errctx_add(errctx, "%s: model objfn must be valid to do profile\n", __func__);

	/* get popmodel and check for errors */
	var popmodel = popmodel_init(source->theta, source->omega, source->sigma, errctx);
	if (errctx->len)
		goto failed;

	/* more error checking or profileconfig */
	var index = args->index - (source->data._offset1 ? 1 : 0);
	switch (args->type) {
		case PROFILE_THETA:
			if (index < 0 || index >= popmodel.ntheta) {
				errctx_add(errctx, "%s: profile theta index (%i) out of bounds", __func__, args->index);
				goto failed;
			}
			break;

		case PROFILE_OMEGA:
			if (index < 0 || index >= popmodel.nomega) {
				errctx_add(errctx, "%s: profile omega index (%i) out of bounds", __func__, args->index);
				goto failed;
			}
			if (popmodel.omegafixed[index][index] == OMEGAFIXED_SAME) {
				errctx_add(errctx, "%s: cannot profile part of omega same block\n", __func__);
				goto failed;
			}
			break;

		case PROFILE_SIGMA:
			if (index < 0 || index >= popmodel.nsigma) {
				errctx_add(errctx, "%s: profile sigma index (%i) out of bounds", __func__, args->index);
				goto failed;
			}
			break;

		default:
			errctx_add(errctx, "%s: profile type (%i) invalid", __func__, args->type);
			goto failed;
	}
	return popmodel;

failed:
	return (POPMODEL) { 0 };
}

static double popmodel_left_value(const POPMODEL* const popmodel,
								 const PROFILECONFIG* const args,
								 const bool _offset1)
{
	var index = args->index - (_offset1 ? 1 : 0);
	switch (args->type) {
		case PROFILE_THETA:
			return popmodel->theta[index];

		case PROFILE_OMEGA:
			return popmodel->omega[index][index];

		case PROFILE_SIGMA:
			return popmodel->sigma[index];

		default:
			return DBL_MAX;
	}
} 

static void popmodel_apply_args(POPMODEL* const popmodel,
								const PROFILECONFIG* const args,
								const bool _offset1)
{
	var index = args->index - (_offset1 ? 1 : 0);
	switch (args->type) {
		case PROFILE_THETA:
			popmodel->theta[index] = args->value;
			popmodel->thetaestim[index] = FIXED;
			break;

		case PROFILE_OMEGA:
			popmodel->omega[index][index] = args->value;
			popmodel->omegafixed[index][index] = OMEGAFIXED_FIXED;
			break;

		case PROFILE_SIGMA:
			popmodel->sigma[index] = args->value;
			popmodel->sigmafixed[index] = 1;
			break;

		default:
			break;
	}
} 

static const char* profile_type(const PROFILECONFIG* const args)
{
	switch (args->type) {
		case PROFILE_THETA:
			return "theta";

		case PROFILE_OMEGA:
			return "omega";

		case PROFILE_SIGMA:
			return "sigma";

		default:
			return "unknown";
	}
}

static void pmx_profile_evaluate_helper(OPENPMX* const ret,
										POPMODEL* const popmodel,
										const OPENPMX* const source,
										PROFILECONFIG* const args,
										const char* name)
{
	/* modify the popmodel */
	popmodel_apply_args(popmodel, args, source->data._offset1);
	memset(&popmodel->result, 0, sizeof(popmodel->result));

	/* reset pmx object and start it from the best point so far */
	pmx_release_state(ret);
	pmx_copy_popmodel(ret, popmodel);
	pmx_ensure_state(ret);
	if (source->state) {
		assert(ret->state->idata.nindivid == source->state->idata.nindivid);
		assert(ret->state->idata.nomega == source->state->idata.nomega);
		idata_set_eta(&ret->state->idata, source->state->idata.individ[0].eta);
	}

	/* figure out the file name, if source has none, we should have none as well */
	char filename[PATH_MAX] = "";
	if (source->filename) {
		let type = profile_type(args);
		snprintf(filename, sizeof(filename), "%s.profile.%s.%i.%s", source->filename, type, args->index, name);
		ret->filename = filename;
	}

	/* do the estimation, with the given filename */
	pmx_estimate(ret, &args->estimate);
	ret->filename = 0; /* dont let filename dangle */
}

static void get_iter_filename(char* filename, const size_t size, const OPENPMX* const source, const PROFILECONFIG* const args, const char* name)
{
	let type = profile_type(args);
	snprintf(filename, size, "%s.profile.%s.%i.%s.iter", source->filename, type, args->index, name);
}

static void write_iter_header(FILE* stream)
{
	fprintf(stream, OPENPMX_SFORMAT OPENPMX_SFORMAT OPENPMX_SFORMAT OPENPMX_SFORMAT OPENPMX_SFORMAT OPENPMX_SFORMAT "\n",
				"type", "index", "value", "objfn", "target", "neval");
}

typedef struct {
	const double left;
	const double right;
	const double objfn_target;
	const OPENPMX* const source;
	PROFILECONFIG* args;
	const char* name;
	POPMODEL* popmodel;
	OPENPMX* ret;
	int neval;
	bool converged;
	FILE* stream;
	bool final_proposal;
	bool test_only;
} PROFILEPARAMS;

static double root_function(double x, void *_params)
{
	var params = (PROFILEPARAMS*)_params;

	/* get the final proposal without evaluation, so we have a proposed
	 * final estimate which makes more sense as an estimate of the
	 * profile crossing point */
	if (params->final_proposal)
		return 0;
	
	let args = params->args;
	let type = profile_type(args);
	let target = params->objfn_target;

	/* construct the new value to test */
	/* sqrt(x) make the root finding more linear and efficient and because
	 * x-range is [0, 1] then there is no issue with domain */
	let delta = params->right - params->left;
	let new_value = params->left + sqrt(x)*delta;
	args->value = new_value;

	double objfn; 
	if (x == 0.) {
		objfn = params->source->result.objfn;
		goto done;
	}

	/* test it and return error */
	var ret = params->ret;
	pmx_profile_evaluate_helper(ret, params->popmodel, params->source, args, params->name);
	objfn = ret->result.objfn;
	params->neval += 1;
	
done:
	if (x != 0. || !params->test_only) {
		if (!params->test_only)
			printf("profile test %s %i %s value %g objfn %.6f target %.6f neval %i\n", 
				   type, args->index, params->name, args->value, objfn, params->objfn_target, params->neval);
		if (params->stream) {
			fprintf(params->stream, OPENPMX_SFORMAT OPENPMX_IFORMAT OPENPMX_FFORMAT 
									OPENPMX_FFORMAT OPENPMX_FFORMAT OPENPMX_IFORMAT "\n",
									type, args->index, args->value, 
									objfn, 
									!params->test_only ? params->objfn_target : 0., 
									params->neval);
		}
	}
   
	let err = objfn - target;
	if (fabs(err) < args->dobjfn_tol) 
		params->converged = true;

	return err;
}

OPENPMX pmx_profile(const OPENPMX* const source, PROFILECONFIG* const args)
{
	ERRCTX errctx = { 0 };
	
/// If .maxeval=1 and .append=true then we go into "test only" mode which 
/// only adds a point to the profile. Then the other fields are ignored.
	let test_only = (args->maxeval == 1 && args->append);

	if (args->dobjfn_tol < 0.)
		fatal(0, "%s: dobjfn_tol (%f) cannot be < 0.\n", __func__, args->dobjfn_tol);
	if (args->maxeval < 0) 
		fatal(0, "%s: maxeval (%i) cannot be < 0.\n", __func__, args->maxeval);

	var popmodel = profile_popmodel_init(source, args, &errctx); 
	if (errctx.len) 
		fatal(0, "%s: %s", __func__, errctx.errmsg);

	let left_value = popmodel_left_value(&popmodel, args, source->data._offset1);
	let right_value = args->value;
	var name = args->name;
	if (!name) {
		name = "right";
		if (right_value < left_value)
			name = "left";
	}
	
	var ret = pmx_copy(source);

	FILE* stream = 0;
	char filename[PATH_MAX] = "";
	let type = profile_type(args);
	if (source->filename) {
		get_iter_filename(filename, sizeof(filename), source, args, name);
		ret.filename = filename;

		/* now open it to add stuff */
		/* try to open to see if it already exists */
		stream = fopen(filename, args->append ? "a" : "w");
		if (!stream)
			fatal(0, "%s: failed to open \"%s\"\n", __func__, filename);
			
		if (!args->append || ftell(stream) == 0)
			write_iter_header(stream);
	}

/// If .dobjfn is 0. then it will be set to the chi-squared distribution 
/// for alpha 0.01 and 1 degree of freedom. This is approximately 6.63.
#define DEFAULT_ALPHA	0.01
#define DEFAULT_DF		1
	if (args->dobjfn == 0.)
		args->dobjfn = gsl_cdf_chisq_Pinv(1.0 - DEFAULT_ALPHA, DEFAULT_DF);
/// The default tolerance is 1 objfn units. This is usually enough that a
/// smoothed line through the profile points is close to the confidence
/// limits.
	if (args->dobjfn_tol <= 0.)
		args->dobjfn_tol = 1.;

/// The default number of root finding evaluations is 10.
	if (args->maxeval == 0)
		args->maxeval = 10;

	/* neval start at one because gsl_root_fsolver_set() will call twice
	 * to bracket the interval. The left side we already take from the 
	 * source so we dont have to reevaluate it, then the right side will
	 * be called */
	var params = (PROFILEPARAMS) {
		.left = left_value,
		.right = right_value,
		.objfn_target = source->result.objfn + args->dobjfn,
		.source = source,
		.args = args,
		.name = name,
		.popmodel = &popmodel,
		.ret = &ret,
		.converged = false,
		.stream = stream,
		.neval = 0,	
		.final_proposal = false,
		.test_only = test_only,
	};
	var F = (gsl_function) {
		.function = root_function,
		.params = &params,
	};

	let vlower = fmin(left_value, args->value);
	let vupper = fmax(left_value, args->value);
	if (!test_only) {
		printf("profile start %s %i lower %g upper %g\n"
			   "profile objfn %.6f target %.6f dobjfn %.6f tol %g\n", 
			   type, args->index, vlower, vupper, 
			   source->result.objfn, params.objfn_target, args->dobjfn, args->dobjfn_tol);
	} else {
		printf("profile test %s %i %s value %g\n", 
			   type, args->index, name, args->value);
	}

	/* bracket in normalised [0, 1] space. The root_function maps back
	 * to [left, right] */
	var s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
	if (!s)
		fatal(0, "%s: profile root solver alloc failed\n", __func__);

	let oldhandler = gsl_set_error_handler_off();
	int status = gsl_root_fsolver_set(s, &F, 0., 1.);
	if (status != GSL_SUCCESS && !test_only) 
		fatal(0, "%s: profile root solver set error %i(%s)\n", __func__, status, gsl_strerror(status));

	/* if we are "test only" then we dont need to iterate the root finder
	 * since the gsl_root_fsolver_set() will have called the bounds and 
	 * thus tested the intended point. */
	if (test_only) {
		printf("profile test complete\n");
		goto cleanup;
	}

	/* iterate the root finder */
	while (1) {
		if (params.neval >= args->maxeval) {
			printf("profile maxeval (%i) reached\n", args->maxeval);
			break;
		}
		status = gsl_root_fsolver_iterate(s);
		if (status != GSL_SUCCESS && status != GSL_CONTINUE) {
			printf("profile root status no continue (%s)\n", gsl_strerror(status));
			break;
		}
		if (params.converged)
			break;
	} 

	/* final estimate, but untested into the output file so the users
	 * gets the interpolated best guess so far. This can only be done
	 * by calling gsl_root_fsolver_iterate() which must call the
	 * root_function(). So we add a flag for this special case and we
	 * can avoid doing a real evaluation */
	if (status == GSL_SUCCESS || status == GSL_CONTINUE) {
		params.final_proposal = true;
		gsl_root_fsolver_iterate(s);
	}

/// If doing profile calculations (not test only) then the .value is 
/// changed to the current best guess as to the point that the profile
/// likelihood crosses the .dobjfn threshold.
	/* reset pmx object and point it to the best guess so far. The
	 * object is emptied so not calling pmx_release_state() is OK. */
	let x = gsl_root_fsolver_root(s);
	let delta = params.right - params.left;
	let new_value = params.left + sqrt(x)*delta;
	args->value = new_value;

	if (params.converged)
		printf("profile converged value %g\n", args->value);
	else
		printf("profile failed value %g\n", args->value);

	if (stream)
		fprintf(stream, OPENPMX_SFORMAT OPENPMX_IFORMAT OPENPMX_FFORMAT 
						OPENPMX_FFORMAT OPENPMX_FFORMAT OPENPMX_IFORMAT "\n",
						type, args->index, args->value, 
						0., params.objfn_target, 0);

	/* cleanup */
cleanup:
	gsl_set_error_handler(oldhandler);
	ret.filename = 0; /* dont let filename dangle */
	pmx_release_state(&ret);

	gsl_root_fsolver_free(s);
	if (stream)
		fclose(stream);
		
	return ret;
}
