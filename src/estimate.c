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
 
/// This file does the outer (stage 2) estimation.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <limits.h>

#include "openpmx.h"
#include "githash.h"
#include "omegainfo.h"
#include "ievaluate.h"
#include "encode.h"
#include "stage1.h"
#include "defines.h"
#include "checkout.h"
#include "print.h"
#include "options.h"
#include "pmxstate.h"
#include "predict.h"
#include "advan/advan.h"
#include "utils/c22.h"
#include "utils/various.h"

#include "buildflags.h"

/*--------------------------------------------------------------------*/
/* different optimizers */
#define OPTIMIZER_OUTER_BOBYQA

typedef struct {
	IDATA* const idata;
	const ADVANFUNCS* const advanfuncs;
	ENCODE test;
	const OPTIONS* const options;
	struct {
		POPMODEL* model;
		double* eta;
	} best;
	struct timespec begin;
	FILE* outstream;
	FILE* extstream;
} STAGE2_PARAMS;

static double get_timestamp(const STAGE2_PARAMS* const params)
{
	struct timespec now;
	clock_gettime(CLOCK_REALTIME, &now);
	return timespec_time_difference(&params->begin, &now) / 1000.;
}

static inline void save_besteta(const STAGE2_PARAMS* const params)
{
	let idata = params->idata;
	let firstindivid = &idata->individ[0];
	memcpy(params->best.eta, firstindivid->eta, idata->nindivid * idata->nomega * sizeof(double));
}

static void print_model_evaluation(const STAGE2_PARAMS* params)
{
	let popmodel = &params->test.popmodel;
	let options = params->options;
	popmodel_eval_information(popmodel,
							  get_timestamp(params),
							  idata_ineval(params->idata, false),
							  options->estimate.details || options->estimate.verbose,
							  params->outstream,
							  params->extstream,
							  0);
}

static void update_best_imodel(const STAGE2_PARAMS* const params)
{
	let popmodel = &params->test.popmodel;
	let options = params->options;
	var best = params->best.model;

	/* update the best estimation and its objective function and function evaluations so far */
	bool print_model = false;
	let dobjfn = popmodel->result.objfn - best->result.objfn;
	if (dobjfn < 0.) { 
/// Each time the objective function improves during estimation the eta
/// values are saved. They will be used as initial values for subsequent
/// optimization.
		*best = *popmodel;
		save_besteta(params);
		print_model = true;
	}

	/* even if we didnt improve, update the number of evaluations of the best model */
	best->result.neval = popmodel->result.neval;

	/* update the user */
	if (options->estimate.verbose)
		print_model = true;
	if (print_model)
		print_model_evaluation(params);
}

static void encode_evaluate(ENCODE* const test,
							IDATA* const idata,
							const ADVANFUNCS* const advanfuncs,
							const OPTIONS* const options)
							
{
	var popmodel = &test->popmodel;
	let omegainfo = &test->omegainfo;
	let nonzero = &omegainfo->nonzero;
	var scatteroptions = (SCATTEROPTIONS) {
		.stage1_order = true,
	};
	scatter_threads(idata, advanfuncs, popmodel, nonzero, options, &scatteroptions, stage1_thread);

	popmodel->result.objfn = idata_objfn(idata, omegainfo->omega_nonzero_lndet);
	popmodel->result.type = OBJFN_CURRENT;
	popmodel->result.neval += 1;
}

static double focei_stage2_evaluate_population_objfn(const long int _xlength,
													 const double _x[static _xlength],
													 void * const data)
{
	/* decode the vector in the right way for the advan and make the popmodel ready to test */
	let params = (STAGE2_PARAMS*) data;
	let idata = params->idata;
	assert(_xlength == params->test.nparam);
	encode_update(&params->test, _x);

	/* do the actual test, this sets the objfn */
	let advanfuncs = params->advanfuncs;
	let options = params->options;
	let popmodel = &params->test.popmodel;
	idata_set_eta(idata, params->best.eta);
	encode_evaluate(&params->test, idata, advanfuncs, options);

	/* update best imodel and inform the user if we improve */
	update_best_imodel(params);

	return popmodel->result.objfn;
}

#ifdef OPTIMIZER_OUTER_BOBYQA
#include "bobyqa/bobyqa.h"
#endif

static bool bobyqa_error(const int retcode, const char* phase, FILE* stream)
{
	if (retcode == BOBYQA_SUCCESS)
		return false;

	static const struct {
		int code;
		const char* msg;
	} errors[] = {
		{ BOBYQA_BAD_NPT,              "NPT is not in the required interval"             },
		{ BOBYQA_TOO_CLOSE,            "insufficient space between the bounds"           },
		{ BOBYQA_ROUNDING_ERRORS,      "too much cancellation in a denominator"          },
		{ BOBYQA_TOO_MANY_EVALUATIONS, "maximum number of function evaluations exceeded" },
		{ BOBYQA_STEP_FAILED,          "a trust region step has failed to reduce Q"      },
	};

	const char* errmsg = "unknown";
	forarray(i, errors) {
		if (errors[i].code == retcode) {
			errmsg = errors[i].msg;
			break;
		}
	}

	warning(stream, "%s BOBYQA error %i: %s\n", phase, retcode, errmsg);
	return true;
}

static bool focei(STAGE2_PARAMS* const params)
{
	assert(params);
	var converged = false;
	let options = params->options;

	let n = params->test.nparam;
	if (n == 0) {
		info(params->outstream, "optim no parameters to estimate\n");
		return true;
	}
	
	var initial = mallocvar(double, n);
	var lower = mallocvar(double, n);
	var upper = mallocvar(double, n);
	forcount(i, n) {
		initial[i] = 0.;
		lower[i] = -1e6; // -DBL_MAX;
		upper[i] = 1e6; // DBL_MAX;
	}
	let maxeval = options->estimate.maxeval;
	let step_initial = options->estimate.step_initial;
	let step_refine = options->estimate.step_refine;
	let step_final = options->estimate.step_final;

#ifdef OPTIMIZER_OUTER_BOBYQA
	let npt = 2*n+1;			/* reccommended */
//	let npt1 = n+2;				/* minimum */
//	let npt2 = (n+1)*(n+2)/2;	/* maximum */
	let iprint = 0;
	let wsize = (npt+5)*(npt+n)+3*n*(n+5)/2 + 10; /* a little bit extra room to be sure */
	let w = mallocvar(double, wsize);

/// Model estimation begins with an intial optimization with large 
/// changes to paramater values with BOBYQA rho values from 
/// step_initial, stopping at step_refine. 
	let best = params->best.model;
	var lastobjfn = best->result.objfn;
	forcount(i, n)
		initial[i] = 0.;
	var neval = maxeval;
	var rhobeg = step_initial;
	var rhoend = step_refine;
	info(params->outstream, "optim rho %g %g\n", rhobeg, rhoend);
	var retcode = bobyqa(n, npt,
						 focei_stage2_evaluate_population_objfn, (void*)params,
						 initial, lower, upper,
						 rhobeg, rhoend,
						 iprint, neval, w);
	bobyqa_error(retcode, "stage2 initial", params->outstream); /* warn error, but try refine anyway */
	neval = maxeval - best->result.neval;

	var timestamp = get_timestamp(params);
	info(params->outstream, "time %.3f neval %i objfn %f\n", timestamp, best->result.neval, best->result.objfn);
	lastobjfn = best->result.objfn;

/// After initial optimization, a smaller search space is used from 
/// rho_refine to rho_final. This is repeated until the change in 
/// objective function between restarts is less than dobjfn.
	if (neval > 1 && step_final < step_refine) {
		while (!converged) {
			if (neval <= 1) {
				timestamp = get_timestamp(params);
				info(params->outstream, "time %.3f optim too few iterations for refine\n", timestamp);
				break;
			}
			encode_offset(&params->test, params->best.model);
			forcount(i, n)
				initial[i] = 0.;

			rhobeg = step_refine;
			rhoend = step_final;
			info(params->outstream, "optim rho %g %g\n", rhobeg, rhoend);
			retcode = bobyqa(n, npt,
							 focei_stage2_evaluate_population_objfn, (void*)params,
							 initial, lower, upper,
							 rhobeg, rhoend,
							 iprint, neval, w);
			neval = maxeval - best->result.neval;
			if (bobyqa_error(retcode, "stage2 refine", params->outstream))
				break;
				
			let dobjfn = best->result.objfn - lastobjfn;
			timestamp = get_timestamp(params);
			info(params->outstream, "time %.3f neval %i objfn %f\n", timestamp, best->result.neval, best->result.objfn);
			lastobjfn = best->result.objfn;

			if (fabs(dobjfn) < fabs(options->estimate.dobjfn))
				converged = true;
		}
	}
	free(w);
#endif

	free(lower);
	free(upper);
	free(initial);

	return converged;
}

static bool stabilize_initial_objfn(STAGE2_PARAMS* params)
{
	let idata = params->idata;
	let advanfuncs = params->advanfuncs;
	let options = params->options;
	let popmodel = &params->test.popmodel;

	info(params->outstream, "stabilize begin\n");
	idata_set_eta(idata, params->best.eta);

	/* very first evaluation */
	encode_evaluate(&params->test, idata, advanfuncs, options);
	var timestamp = get_timestamp(params);
	info(params->outstream,
		 "time %.3f neval %i objfn %f\n",
		 timestamp,
		 popmodel->result.neval,
		 popmodel->result.objfn);
	if (!isfinite(popmodel->result.objfn))
		fatal(params->outstream, "objective function not finite\n");
	
/// At initial objective function calculation the results may not be 
/// stable with recalculating resulting in a different objective 
/// function value. This seems likely due to the changing of the 
/// initial conditions for each Stage 1 due to updating of the best eta 
/// value. So for the first evaluation, the objective function value is
/// recalculated until it is stabel to 0.01 or 10 iterations. If this
/// stabilization phase fails or requires many iterations then there may
/// be numerical stability issues with the model and data.
	int i;
	let maxiter = 10;
	for (i=0; i<maxiter; i++) {
		let lastobjfn = popmodel->result.objfn;
		encode_evaluate(&params->test, idata, advanfuncs, options);
		print_model_evaluation(params);

		/* are we stable? */
		let dobjfn = popmodel->result.objfn - lastobjfn;
		if (fabs(dobjfn) < 0.01)
			break;
	}

	/* warn if we are not stable after 10 iterations */
	if (i >= maxiter) 
		warning(params->outstream, "stabilize failed\n");

	/*save the popmodel shouldnt actually have been changed */
	save_besteta(params);
	*params->best.model = *popmodel;

	return (i < maxiter);
}

static void focei_popmodel_stage2(STAGE2_PARAMS* params)
{
	let idata = params->idata;
	let options = params->options;
	let outstream = params->outstream;

	/* cleanup previous runs */
	/* at each estimate or evaluate, the pred and state gets reset to zero */
	let ndata = idata->ndata;
	var firstindivid = &idata->individ[0];
	memset(firstindivid->pred, 0, ndata * sizeof(double));
	memset(firstindivid->istate, 0, ndata * idata->nstate * sizeof(double));

	idata_free_simerr(idata);
	if (params->options->estimate.stage1.icov_resample)
		idata_alloc_icovresample(idata);
	else
		idata_free_icovresample(idata);

	/* first run so we can set objective function and yhat. */
	/* this is the first evaluation, the encode_offset at the best model is
	 * already done. We may have to evaluate several times for a stable objfn */
	/* after the iterations, we save the besteta which we use to start with later */
	if (!stabilize_initial_objfn(params))
		warning(outstream, "initial objfn not stable\n");

	/* make sure best has the objfn of what we just stabilized */
	var best = params->best.model;
	*best = params->test.popmodel;
	encode_offset(&params->test, best); /* this probably does nothing */

	/* call the underlying advan to optimize and find the best imodel */
	let maxeval = options->estimate.maxeval;
	best->result.type = OBJFN_EVALUATE;
	if (maxeval > 1) {
		let fullconverge = focei(params); /* TODO: this should be renamed I think */
		if (fullconverge) 
			best->result.type = OBJFN_FINAL;

		/* at end we encode the best so far */
		encode_offset(&params->test, best);
	}
}

typedef enum {
	OUTFILE_HEADER_EVALUATE,
	OUTFILE_HEADER_ESTIMATE,
	OUTFILE_HEADER_RESAMPLE
} OUTFILE_TYPE;

static void outfile_header(FILE* f2,
						   const ADVANFUNCS* const advanfuncs,
						   const IDATA* const idata,
						   const OPTIONS* const options,
						   const OUTFILE_TYPE outfile_type)
{
	assert(advanfuncs);
	assert(idata);
	assert(options);
	
	info(f2, "OpenPMX %i.%i.%i hash %s\n", OPENPMX_VERSION_MAJOR, OPENPMX_VERSION_MINOR, OPENPMX_VERSION_RELEASE, OPENPMX_GITHASH);

#if defined(OPENPMX_PARALLEL_PTHREADS)
	char parallel_message[] = "pthread";
#elif defined(OPENPMX_PARALLEL_OPENMP)
	char parallel_message[] = "openmp";
#elif defined(OPENPMX_PARALLEL_SINGLETHREAD)
	char parallel_message[] = "singlethread";
#else
#error no parallel processing defined
#endif

#if defined(OPENPMX_SERVER)
	char server_message[] = "server";
#elif defined(OPENPMX_NOTSERVER)
	char server_message[] = "notserver";
#else
#error no server option defined
#endif

	info(f2, "config %s %i %s %s\n", 
		parallel_message, options->nthread, server_message, OPENPMX_INSTALL_PLATFORM);

	info(f2, "data records %i observations %i individuals %i\n", advanfuncs->recordinfo.dataconfig->nrecords, idata->nobs, idata->nindivid);
	info(f2, "data table offset %s\n", advanfuncs->recordinfo.dataconfig->_offset1 ? "true" : "false");

	advanfuncs->info(advanfuncs, stdout);
	if (f2)
		advanfuncs->info(advanfuncs, f2);

	if (outfile_type == OUTFILE_HEADER_RESAMPLE)
		info(f2, "resample seed %i\n", options->simulate.seed);
}

static FILE* estimate_results_fopen(const char* name,
									const char* ext)
{
	char fname[PATH_MAX];
	snprintf(fname, sizeof(fname), "%s%s", name, ext);
	return fopen(fname, "w");
}

static STAGE2_PARAMS stage2_params_init(const char* filename,
										IDATA* const idata,
										const ADVANFUNCS* const advanfuncs,
										POPMODEL* const popmodel,
										const OPTIONS* const options)
{
	/* setup output logging */
	FILE* outstream = 0;
	FILE* extstream = 0;
	if (filename) {
		outstream = estimate_results_fopen(filename, OPENPMX_OUTFILE);
		if (!outstream)
			fatal(0, "%s: could not open file \"%s%s\"\n", __func__, filename, OPENPMX_OUTFILE);
		extstream = estimate_results_fopen(filename, OPENPMX_EXTFILE);
		if (!extstream)
			fatal(outstream, "%s: could not open file %s extension %s\n", __func__, filename, OPENPMX_EXTFILE);
	}

	var params = (STAGE2_PARAMS) {
		.idata = idata,
		.advanfuncs = advanfuncs,
		.test = encode_init(popmodel),	/* makes a copy of popmodel to work with */
		.options = options,
		.best = {
			.model = popmodel,			/* best so far will be saved at the caller */
			.eta = callocvar(double, idata->nindivid * idata->nomega),					
		},
		.begin = { }, 					/* set after initialization */
		.outstream = outstream,
		.extstream = extstream,
	};
	clock_gettime(CLOCK_REALTIME, &params.begin);
	if (idata->nindivid <= 0)
		fatal(outstream, "%s: optim cannot estimate, no individuals\n", __func__);

	/* save the current, but with objfn to invalid and nparam as well */
	*params.best.model = params.test.popmodel;

	/* make sure 0 values reproduce the current model */
	encode_offset(&params.test, params.best.model);

	return params;
}

static void stage2_params_cleanup(STAGE2_PARAMS *params)
{
	if (params->outstream)
		fclose(params->outstream);
	if (params->extstream)
		fclose(params->extstream);
	free(params->best.eta);
}

static void estimate_popmodel(const char* filename,
							  IDATA* const idata,
							  const ADVANFUNCS* const advanfuncs,
							  POPMODEL* const popmodel,
							  const OPTIONS* const options)
{
	const char* message = "Evaluate only";
	var outfile_type = OUTFILE_HEADER_EVALUATE;
	let maxeval = options->estimate.maxeval;
	popmodel->result.type = OBJFN_CURRENT;
	if (maxeval > 1) {
#ifdef OPTIMIZER_OUTER_BOBYQA
		message = "FOCEI BOBYQA";
#endif
		outfile_type = OUTFILE_HEADER_ESTIMATE;
	}

	/* setup the minimizer */
	var params = stage2_params_init(filename, idata, advanfuncs, popmodel, options);

/// At start of estimation the header of the ext file is written.
	/* do some logging */
	/* start extfile and other headers, rest will be saved during iterations */
	outfile_header(params.outstream, advanfuncs, idata, options, outfile_type);
	var _offset1 = advanfuncs->recordinfo.dataconfig->_offset1;
	if (params.extstream) 
		extfile_header(params.extstream, params.best.model, _offset1);
	if (maxeval > 1)
		info(params.outstream, "optim %s\n", message);

	/* update results to screen and log */
	popmodel_information(params.outstream, popmodel, -DBL_MAX);

/// Before estimation a data checkout is done (see checkout.c and 
/// ievaluate.c) to detect various errors. 
	idata_ineval(idata, true);
	idata_checkout(idata, advanfuncs, popmodel, options, params.outstream);

	/* actually do the stage 2 estimation */
	/* best.model in params saves the best at the caller so we dont
	 * need to save it ourselves */
	focei_popmodel_stage2(&params);

	/* update results to screen and log */
	let timestamp = get_timestamp(&params);
	popmodel_information(params.outstream,
						 popmodel,
						 timestamp);

/// At the end of estimation the phi file is written and a trailer is 
/// put onto the ext file. Also the yhat file is written with prediction
/// varaiables including pred, yhat, yhatvar, and the state.
	/* update results in tables */
	if (filename) {
		table_phi_idata(filename, idata, _offset1);

		/* TODO: could this be optional, it might be burdensome for 
		 * large models because it forces another advancing throught the 
		 * records */
		idata_predict_pred(idata, advanfuncs, popmodel, options);
		table_yhat_idata(filename, idata, _offset1);

		if (params.extstream) {
			let ineval = idata_ineval(idata, false);
			extfile_trailer(params.extstream, popmodel, timestamp, ineval);
		}
		if (options->estimate.stage1.icov_resample)
			table_icov_resample_idata(filename, idata, _offset1);
	}
	
	/* cleanup */
	stage2_params_cleanup(&params);
}

void pmx_estimate(OPENPMX* pmx, ESTIMCONFIG* const estimate)
{
	pmxstate_ensure(pmx);
	var pstate = pmx->state;

	var options = options_init(pmx);
	if (estimate) {
		options.estimate = estimconfig_default(estimate);
		*estimate = options.estimate;
	}

	ERRCTX errctx = { 0 };
	var popmodel = popmodel_init(pmx->theta, pmx->omega, pmx->sigma, &errctx);
	if (errctx.len)
		fatal(0, "%s: %s", __func__, errctx.errmsg);

	estimate_popmodel(pmx->filename,
					  &pstate->idata,
					  pstate->advanfuncs,
					  &popmodel,
					  &options);

	pmx_popmodel_writeback(pmx, &popmodel);
}

/// Evaluation is the same as estimation but with maxeval=0.

void pmx_evaluate(OPENPMX* pmx, STAGE1CONFIG* const stage1)
{
	pmxstate_ensure(pmx);
	var pstate = pmx->state;

	var options = options_init(pmx);
	if (stage1) {
		options.estimate.stage1 = stage1config_default(stage1);
		*stage1 = options.estimate.stage1;
	}
	options.estimate.maxeval = 1;

	ERRCTX errctx = { 0 };
	var popmodel = popmodel_init(pmx->theta, pmx->omega, pmx->sigma, &errctx);
	if (errctx.len)
		fatal(0, "%s: %s", __func__, errctx.errmsg);

	estimate_popmodel(pmx->filename,
					  &pstate->idata,
					  pstate->advanfuncs,
					  &popmodel,
					  &options);

	pmx_popmodel_writeback(pmx, &popmodel);
}

void pmx_fastestimate(OPENPMX* pmx, ESTIMCONFIG* const estimate)
{
	pmxstate_ensure(pmx);
	var pstate = pmx->state;

	var options = options_init(pmx);
	if (estimate) {
		options.estimate = estimconfig_default(estimate);
		*estimate = options.estimate;
	}
	options.estimate.step_initial = 1.;
	options.estimate.step_refine = 0.1;
	options.estimate.step_final = 0.01;

	ERRCTX errctx = { 0 };
	var popmodel = popmodel_init(pmx->theta, pmx->omega, pmx->sigma, &errctx);
	if (errctx.len)
		fatal(0, "%s: %s", __func__, errctx.errmsg);

	estimate_popmodel(pmx->filename,
					  &pstate->idata,
					  pstate->advanfuncs,
					  &popmodel,
					  &options);

	pmx_popmodel_writeback(pmx, &popmodel);
}


