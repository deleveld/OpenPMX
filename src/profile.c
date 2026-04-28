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

static POPMODEL pmx_profile_initial_popmodel(const OPENPMX* const source,
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
			if (index < 0 || index >= popmodel.ntheta)
				errctx_add(errctx, "%s: profile theta index (%i) out of bounds", __func__, args->index);
			break;

		case PROFILE_OMEGA:
			if (index < 0 || index >= popmodel.nomega)
				errctx_add(errctx, "%s: profile omega index (%i) out of bounds", __func__, args->index);
			if (popmodel.omegafixed[index][index] == OMEGAFIXED_SAME)
				errctx_add(errctx, "%s: cannot profile part of omega same block\n", __func__);
			break;

		case PROFILE_SIGMA:
			if (index < 0 || index >= popmodel.nsigma)
				errctx_add(errctx, "%s: profile sigma index (%i) out of bounds", __func__, args->index);
			break;

		default:
			errctx_add(errctx, "%s: profile type (%i) invalid", __func__, args->type);
	}
	return popmodel;

failed:
	return (POPMODEL) { 0 };
}

static const char profile_theta[] = "theta";
static const char profile_omega[] = "omega";
static const char profile_sigma[] = "sigma";

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
			return profile_theta;

		case PROFILE_OMEGA:
			return profile_omega;

		case PROFILE_SIGMA:
			return profile_sigma;

		default:
			return 0;
	}
}

static void pmx_profile_evaluate_helper(OPENPMX* const ret,
										POPMODEL* const popmodel,
										const OPENPMX* const source,
										PROFILECONFIG* const args,
										const char* direction)
{
	/* modify the popmodel */
	popmodel_apply_args(popmodel, args, source->data._offset1);

	/* reset pmx object and start it from the best point so far */
	pmx_cleanup(ret);
	pmx_popmodel_writeback(ret, popmodel);
	pmxstate_ensure(ret);
	if (source->state) {
		assert(ret->state->idata.nindivid == source->state->idata.nindivid);
		assert(ret->state->idata.nomega == source->state->idata.nomega);
		idata_set_eta(&ret->state->idata, source->state->idata.individ[0].eta);
	}

	/* figure out the file name, if source has none, we should have none as well */
	char filename[PATH_MAX] = "";
	if (source->filename) {
		let type = profile_type(args);
		snprintf(filename, sizeof(filename), "%s.profile.%s.%i.%s", source->filename, type, args->index, direction);
		ret->filename = filename;
	}

	/* do the estimation */
	pmx_estimate(ret, &args->estimate);
	ret->filename = 0; /* dont let filename dangle */

	/* pass results back by modifying config */
	args->dobjfn = ret->result.objfn - source->result.objfn;
}

OPENPMX pmx_profile_evaluate(const OPENPMX* const source, PROFILECONFIG* const args)
{
	ERRCTX errctx = { 0 };

	var popmodel = pmx_profile_initial_popmodel(source, args, &errctx); 
	if (errctx.len) {
		warning(0, "%s: %s", __func__, errctx.errmsg);
		return (OPENPMX) { 0 };
	}

	let left_value = popmodel_left_value(&popmodel, args, source->data._offset1);
	var direction = "right";
	if (args->value < left_value)
		direction = "left";

	var ret = pmx_copy(source);
	pmx_profile_evaluate_helper(&ret, &popmodel, source, args, direction);
	return ret;
}

typedef struct {
	const double left;
	const double right;
	const double dobjfn_target;
	const double objfn_target;
	const OPENPMX* const source;
	PROFILECONFIG* args;
	const char* direction;
	POPMODEL* popmodel;
	OPENPMX* ret;
	int neval;
	bool converged;
	FILE* stream;
} PROFILEPARAMS;

static double root_function(double x, void *_params)
{
	var params = (PROFILEPARAMS*)_params;
	let args = params->args;

	if (x == 0.) 
		return -params->dobjfn_target;
	
	/* construct the new value to test */
	/* sqrt(x) make the root finding more linear and efficient and because
	 * x-range is [0 1] then there is no issue with domain */
	let delta = params->right - params->left;
	let new_value = params->left + sqrt(x)*delta;
	args->value = new_value;

	/* test it and return error */
	var ret = params->ret;
	pmx_profile_evaluate_helper(ret, params->popmodel, params->source, args, params->direction);
	params->neval += 1;

	let err = ret->result.objfn - params->objfn_target;
	info(params->stream, "profile value %f objfn %.6f error %.6f\n",
						 new_value, ret->result.objfn, err);
   
	if (fabs(err) < args->dobjfn_tol) 
		params->converged = true;

	return err;
}

OPENPMX pmx_profile(const OPENPMX* const source, PROFILECONFIG* const args)
{
	ERRCTX errctx = { 0 };

	var popmodel = pmx_profile_initial_popmodel(source, args, &errctx); 
	if (errctx.len) {
		warning(0, "%s: %s", __func__, errctx.errmsg);
		return (OPENPMX) { 0 };
	}

	let left_value = popmodel_left_value(&popmodel, args, source->data._offset1);
	var direction = "right";
	if (args->value < left_value)
		direction = "left";

	var ret = pmx_copy(source);

	FILE* stream = 0;
	char filename[PATH_MAX] = "";
	let type = profile_type(args);
	if (source->filename) {
		snprintf(filename, sizeof(filename), "%s.profile.%s.%i.%s.iter", source->filename, type, args->index, direction);
		ret.filename = filename;

		stream = fopen(filename, "w");
		if (!stream)
			fatal(0, "%s: failed to open \"%s\"\n", __func__, filename);
			
		fprintf(stream, "OpenPMX %i.%i.%i hash %s\n", OPENPMX_VERSION_MAJOR, OPENPMX_VERSION_MINOR, OPENPMX_VERSION_RELEASE, OPENPMX_GITHASH);
	}

	/* default tolerance is 0.01 objfn units */
	let default_alpha = 0.01;
	let default_df = 1.;
	if (args->dobjfn <= 0.)
		args->dobjfn = gsl_cdf_chisq_Pinv(1.0 - default_alpha, default_df);
	if (args->dobjfn_tol <= 0.)
		args->dobjfn_tol = 0.1;
	if (args->maxeval <= 0)
		args->maxeval = 10;

	var params = (PROFILEPARAMS) {
		.left = left_value,
		.right = args->value,
		.dobjfn_target = args->dobjfn,
		.objfn_target = source->result.objfn + args->dobjfn,
		.source = source,
		.args = args,
		.direction = direction,
		.popmodel = &popmodel,
		.ret = &ret,
		.converged = false,
		.stream = stream,
		.neval = 0,
	};
	var F = (gsl_function) {
		.function = root_function,
		.params = &params,
	};

	let vlower = fmin(left_value, args->value);
	let vupper = fmax(left_value, args->value);
	info(stream, "profile start %s %i lower %g upper %g\nprofile target %.6f dobjfn %.6f tol %g\n", type, args->index, vlower, vupper, params.objfn_target, args->dobjfn, args->dobjfn_tol);

	/* bracket in normalised [0, 1] space. The root_function maps back
	 * to [left, right] */
	var s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
	if (!s)
		fatal(0, "%s: root solver alloc failed\n", __func__);

	int status = gsl_root_fsolver_set(s, &F, 0., 1.);
	if (status != GSL_SUCCESS)
		fatal(0, "%s: root solver bracket invalid (%s)", __func__,
			  gsl_strerror(status));

	do {
		status = gsl_root_fsolver_iterate(s);
		if (status != GSL_SUCCESS && status != GSL_CONTINUE) {
			info(stream, "profile root status no continue (%s)\n", gsl_strerror(status));
			break;
		}
		if (params.neval >= args->maxeval) {
			info(stream, "profile maxeval (%i) reached\n", args->maxeval);
			break;
		}
	} while (!params.converged);

	if (params.converged)
		info(stream, "profile converged dobjfn %.6f neval %i\n", args->dobjfn, params.neval);
	else
		info(stream, "profile failed dobjfn %.6f neval %i\n", args->dobjfn, params.neval);

	/* cleanup */
	gsl_root_fsolver_free(s);

	ret.filename = 0; /* dont let filename dangle */
	if (stream)
		fclose(stream);

	return ret;
}
