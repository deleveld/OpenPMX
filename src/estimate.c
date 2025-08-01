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
#include "advan/advan.h"
#include "utils/c22.h"
#include "utils/various.h"
#include "openpmx_internal.h"

#include "openpmx_compile_options.h"

/*--------------------------------------------------------------------*/
/* different optimizers */

#define OPTIMIZER_OUTER_BOBYQA
//#define OPTIMIZER_OUTER_GSL_NELDERMEAD
//#define OPTIMIZER_OUTER_GSL_BFGS

static void pmx_update_from_popmodel(OPENPMX* const pmx, const POPMODEL* const popmodel)
{
	let theta = popmodel->theta;
	let ntheta = popmodel->ntheta;
	let thetaestim = popmodel->thetaestim;
	forcount(i, ntheta) {
		if (thetaestim[i] == ESTIMATE)
			pmx->theta[i].value = theta[i];
	}

	let omega = popmodel->omega;
	let omegafixed = popmodel->omegafixed;
	var offset = 0;
	forcount(k, popmodel->nblock) {
		let ndim = popmodel->blockdim[k];
		let type = popmodel->blocktype[k];

		assert((int)pmx->omega[k].type == type);
		assert(pmx->omega[k].ndim == ndim);

		switch (type) {

			case OMEGA_DIAG: {
				forcount(i, ndim) {
					double v = omega[offset + i][offset + i];
					let f = omegafixed[offset + i][offset + i];
					if (f != 0 && v != 0.)
						v = -fabs(v);
					pmx->omega[k].values[i] = v;
				}
				break;
			}
			case OMEGA_BLOCK: {
				int blocki = 0;
				forcount(i, ndim) {
					forcount(j, i + 1) {
						double v = omega[offset + i][offset + j];
						let f = omegafixed[offset + i][offset + j];
						if (f != 0 && v != 0.)
							v = -fabs(v);
						pmx->omega[k].values[blocki] = v;
						blocki++;
					}
				}
				break;
			}
			case OMEGA_SAME: {
				/* we dont need to do anything here since the block is already */
				break;
			}
			default:
				fatal(0, "Invalid block type (%i) if size %i\n", type, ndim);
				break;
		}
		offset += ndim;
	}

	let sigma = popmodel->sigma;
	let nsigma = popmodel->nsigma;
	let sigmafixed = popmodel->sigmafixed;
	forcount(i, nsigma) {
		var v = sigma[i];
		var f = sigmafixed[i];
		if (f == 0)
			pmx->sigma[i] = v;
	}

	pmx->result = popmodel->result;
}

static double objfn(const IDATA* const idata,
				    const OMEGAINFO* const omegainfo)
{
	let nindivid = idata->nindivid;

	/* doing the sum by type makes sure we are adding numbers of comparable
	 * magnitude which helps with accuracy */
	double objfn1 = 0.;
	double objfn2 = 0.;
	double objfn3 = 0.;
	double objfn4 = 0.;
	double objfn5 = 0.;
	forcount(k, nindivid) {
		let individ = &idata->individ[k];

		let term1 = individ->obs_lndet;
		let term2 = individ->obs_min2ll;
		let term3 = individ->eta_min2ll;
		var term4 = omegainfo->omega_nonzero_lndet;
 		let term5 = individ->icov_lndet;
		if (individ->nobs == 0) {
			assert(term1 == 0.);
			assert(term2 == 0.);
			assert(term3 == 0.);
			term4 = 0.;				/* population term does not count if no observations in the individual */
			assert(term5 == 0.);
		}

		assert(isfinite(term1) == 1);
		assert(isfinite(term2) == 1);
		assert(isfinite(term3) == 1);
		assert(isfinite(term4) == 1);
		assert(isfinite(term5) == 1);

		let iobjfn = term1 + term2 + term3 + term4 + term5;
		individ->iobjfn = iobjfn;

		objfn1 += term1;
		objfn2 += term2;
		objfn3 += term3;
		objfn4 += term4;
		objfn5 += term5;
	}
	let objfn = objfn1 + objfn2 + objfn3 + objfn4 + objfn5;
	assert(isfinite(objfn) == 1);
	return objfn;
}

typedef struct {
	IDATA* const idata;
	const ADVANFUNCS* const advanfuncs;
	ENCODE test;
	const OPTIONS* const options;

	int neval;
	POPMODEL best;
	double* besteta;
	const int neta;

	struct timespec begin;
	FILE* outstream;
	FILE* extstream;
	const char* filename;
} STAGE2_PARAMS;

static double get_timestamp(const STAGE2_PARAMS* const params)
{
	struct timespec now;
	clock_gettime(CLOCK_REALTIME, &now);
	return timespec_time_difference(&params->begin, &now) / 1000.;
}

static void update_best_imodel(const STAGE2_PARAMS* const params,
							   POPMODEL* const best)
{
	let idata = params->idata;
	let popmodel = &params->test.popmodel;
	let options = params->options;

	const POPMODEL* improved_model = 0;
	let dobjfn = popmodel->result.objfn - best->result.objfn;
	if (dobjfn < 0.) { 
		/* update the best estimation and its objective function and function evaluations so far */
		*best = *popmodel;
		improved_model = best;

		/* save the best eta so we can keep restarting there for speed and hopefully some consistancy */
		let firstindivid = &idata->individ[0];
		memcpy(params->besteta, firstindivid->eta, params->neta * sizeof(double));
	}

	/* update the user */
	if (options->estimate.verbose)
		improved_model = popmodel;
	if (improved_model) {
		let runtime_s = get_timestamp(params);
		let ineval = idata_ineval(idata, false);
		var outstream = (options->estimate.progress) ? (params->outstream) : 0;
		var extstream = (options->estimate.progress) ? (params->extstream) : 0;
		popmodel_eval_information(improved_model,
								  runtime_s,
								  ineval,
								  options->estimate.details || options->estimate.verbose,
								  outstream,
								  extstream,
								  0);
	}
}

static void reset_eta(IDATA* const idata, const double* eta)
{
	assert(eta);
	var firstindivid = &idata->individ[0];
	memcpy(firstindivid->eta, eta, idata->nindivid * idata->nomega * sizeof(double));
}

static void encode_evaluate(ENCODE* const test,
							IDATA* const idata,
							const ADVANFUNCS* const advanfuncs,
							const OPTIONS* const options)
							
{
	var popmodel = &test->popmodel;
	let omegainfo = &test->omegainfo;
	let nonzero = &omegainfo->nonzero;
	SCATTEROPTIONS scatteroptions = { };
	scatteroptions.stage1_order = true;
	scatter_threads(idata, advanfuncs, popmodel, nonzero, options, &scatteroptions, stage1_thread);

	popmodel->result = (PMXRESULT) {
		.objfn = objfn(idata, omegainfo),
		.type = OBJFN_CURRENT,
		.neval = 0 };
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
	reset_eta(idata, params->besteta);
	encode_evaluate(&params->test, idata, advanfuncs, options);
	params->neval += 1;
	popmodel->result.neval = params->neval;

	/* update best imodel and inform the user if we improve */
	update_best_imodel(params, &params->best);

	return popmodel->result.objfn;
}

#ifdef OPTIMIZER_OUTER_BOBYQA
#include "bobyqa/bobyqa.h"
#endif

#ifdef OPTIMIZER_OUTER_GSL_NELDERMEAD
#include <gsl/gsl_multimin.h>
static double gsl_stage2_objfn(const gsl_vector* v, void *params)
{
	return focei_stage2_evaluate_population_objfn(v->size, v->data, params);
}
#endif

#ifdef OPTIMIZER_OUTER_GSL_BFGS
#include <gsl/gsl_multimin.h>
static double foce_stage2_f(const gsl_vector* const v, void* const params)
{
	return focei_stage2_evaluate_population_objfn(v->size, v->data, params);
}

static void foce_stage2_df(const gsl_vector * const x, void* const _params, gsl_vector* const J)
{
	let params = (const STAGE2_PARAMS *) _params;
	let step_size = params->options->estimate.step_initial;
	let n = params->test.nparam;
	var xx = gsl_vector_alloc(n);

	let h = step_size;
	forcount(i, n) {
		let v = gsl_vector_get(x, i);
		gsl_vector_memcpy(xx, x);

		let above1 = v + h;
		gsl_vector_set(xx, i, above1);
		let f_above1 = foce_stage2_f(xx, _params);

		let below1 = v - h;
		gsl_vector_set(xx, i, below1);
		let f_below1 = foce_stage2_f(xx, _params);

		let deriv = (f_above1 - f_below1) / (above1 - below1);	/* arbitrary method */

		gsl_vector_set(J, i, deriv);
	}
	gsl_vector_free(xx);
}
static void foce_stage2_fdf(const gsl_vector* const x,
							void* data,
							double* const f,
							gsl_vector* const J)
{
	*f = foce_stage2_f(x, data);
	foce_stage2_df(x, data, J);
}
#endif

static const char* focei(STAGE2_PARAMS* const params)
{
	if (!params)
		return "FOCEI BOBYQA";

	let options = params->options;

	let n = params->test.nparam;
	var initial = mallocvar(double, n);
	var upper = mallocvar(double, n);
	var lower = mallocvar(double, n);
	forcount(i, n) {
		initial[i] = 0.;
		lower[i] = -DBL_MAX;
		upper[i] = DBL_MAX;
	}
	let maxeval = options->estimate.maxeval;
	let step_initial = options->estimate.step_initial;
	let step_refine = options->estimate.step_refine;
	let step_final = options->estimate.step_final;

#ifdef OPTIMIZER_OUTER_GSL_NELDERMEAD
	/* Initialize method and iterate */
	var minex_func = (gsl_multimin_function) {
		.n = n,
		.f = gsl_stage2_objfn,
		.params = params,
	};

	var gsl_initial = gsl_vector_view_array(initial, n);
	let ss = gsl_vector_alloc(n);
	gsl_vector_set_all(ss, step_initial); /* Starting step size */

	var s = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, n);
	gsl_multimin_fminimizer_set(s, &minex_func, &gsl_initial.vector, ss);
	int iter = 0;
	int status;
	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status)
			break;

		let size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, step_final);
	}
	while (status == GSL_CONTINUE && iter < maxeval);

	gsl_multimin_fminimizer_free(s);
	gsl_vector_free(ss);
#endif

#ifdef OPTIMIZER_OUTER_GSL_BFGS
	var my_func = (gsl_multimin_function_fdf) {
		.f = foce_stage2_f,
		.df = foce_stage2_df, /* TODO: consider http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/ */
		.fdf = foce_stage2_fdf,
		.n = (size_t)n,
		.params = (void*)params,
	};

	/* setup and run minimizer */
	/* arbitrary values */
	var prevobjfn = DBL_MAX;
	while (1) {
		var gsl_initial = gsl_vector_view_array(initial, n);
		var s = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2, n);
		gsl_multimin_fdfminimizer_set(s, &my_func, &gsl_initial.vector, step_refine, 0.1);

		/* multimin iterate */
		int status;
		while (1) {
			status = gsl_multimin_fdfminimizer_iterate(s);
			if (status == GSL_ENOPROG)
				break;
			let grad = gsl_multimin_fdfminimizer_gradient(s);
			status = gsl_multimin_test_gradient(grad, step_final);
			if (status !=  GSL_CONTINUE)
				break;
			if (params->neval >=  maxeval) 
				break;
		}
		gsl_multimin_fdfminimizer_free(s);

		let objfn = params->best.result.objfn;
		if (prevobjfn - objfn < fabs(options->estimate.dobjfn))
			break;
		prevobjfn = objfn;
	}
#endif

#ifdef OPTIMIZER_OUTER_BOBYQA
	let npt = 2*n+1;			/* reccommended */
//	let npt1 = n+2;				/* minimum */
//	let npt2 = (n+1)*(n+2)/2;	/* maximum */
	let iprint = 0;
	let wsize = (npt+5)*(npt+n)+3*n*(n+5)/2 + 10; /* a little bit extra room to be sure */
	let w = mallocvar(double, wsize);

	let best = &params->best;
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
	if (retcode != BOBYQA_SUCCESS)
		info(params->outstream, "initial BOBYQA error code %i\n", retcode);

	var timestamp = get_timestamp(params);
	info(params->outstream, "time %.3f neval %i objfn %f\n", timestamp, params->neval, best->result.objfn);
	lastobjfn = best->result.objfn;

	if (step_final < step_refine) {
		var dobjfn = best->result.objfn - lastobjfn;
		do {
			encode_offset(&params->test, &params->best);
			forcount(i, n)
				initial[i] = 0.;

			neval = maxeval - params->neval;
			rhobeg = step_refine;
			rhoend = step_final;
			info(params->outstream, "optim rho %g %g\n", rhobeg, rhoend);
			retcode = bobyqa(n, npt,
							 focei_stage2_evaluate_population_objfn, (void*)params,
							 initial, lower, upper,
							 rhobeg, rhoend,
							 iprint, neval, w);
			if (retcode != BOBYQA_SUCCESS)
				info(params->outstream, "refine BOBYQA error code %i\n", retcode);

			dobjfn = best->result.objfn - lastobjfn;
			var timestamp = get_timestamp(params);
			info(params->outstream, "time %.3f neval %i objfn %f\n", timestamp, params->neval, best->result.objfn);
			lastobjfn = best->result.objfn;
		}
		while (dobjfn < -1.*fabs(options->estimate.dobjfn));
	}
	info(params->outstream, "optim rho %g\n", rhoend);

	free(w);
#endif

	free(lower);
	free(upper);
	free(initial);

	return 0;
}

static void print_model(STAGE2_PARAMS* params)
{
	let idata = params->idata;
	let popmodel = &params->test.popmodel;
	let options = params->options;

//	char message[128] = "";
//	if (popmodel->result.objfn != DBL_MAX && params->best.result.objfn != DBL_MAX)
//		sprintf(message, " objfn %f", popmodel->result.objfn);

	let runtime_s = get_timestamp(params);
	let ineval = idata_ineval(idata, false);
	var outstream = (options->estimate.progress) ? (params->outstream) : 0;
	var extstream = (options->estimate.progress) ? (params->extstream) : 0;
	popmodel_eval_information(popmodel,
							  runtime_s,
							  ineval,
							  options->estimate.details || options->estimate.verbose,
							  outstream,
							  extstream,
							  0);
}

static bool stabilize_model(STAGE2_PARAMS* params)
{
	let idata = params->idata;
	let advanfuncs = params->advanfuncs;
	let options = params->options;
	let popmodel = &params->test.popmodel;

	info(params->outstream, "stabilize begin\n");
	reset_eta(idata, params->besteta);

	/* very first evaluation */
	encode_evaluate(&params->test, idata, advanfuncs, options);
	params->neval += 1;
	var timestamp = get_timestamp(params);
	info(params->outstream, "time %.3f neval %i objfn %f\n", timestamp, params->neval, popmodel->result.objfn);

	var done = false;
	var niter = 0;
	while (!done) {
		let lastobjfn = popmodel->result.objfn;
		encode_evaluate(&params->test, idata, advanfuncs, options);
		params->neval += 1;
		++niter;
		popmodel->result.neval = params->neval;

		print_model(params);

		/* are we stable? */
		/* besteta we update later, individual eta values keep getting updated */
		/* warn if we are not stable after 10 iterations */
		let dobjfn = popmodel->result.objfn - lastobjfn;
		if (fabs(dobjfn) < 0.01) 
			done = true;
		else if (niter >= 10)
			break;
			
		params->best = *popmodel;
	}
	info(params->outstream, "stabilize iter %i\n", niter);
	timestamp = get_timestamp(params);
	info(params->outstream, "time %.3f neval %i objfn %f\n", timestamp, params->neval, popmodel->result.objfn);

	var firstindivid = &idata->individ[0];
	memcpy(params->besteta, firstindivid->eta, idata->nindivid * idata->nomega * sizeof(double));

	return done;
}

#if 0
#include <gsl/gsl_spline.h>
static void evaluate_gradient(STAGE2_PARAMS* params)
{
	/* http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/#noiserobust_2 */
	
	let idata = params->idata;
	let advanfuncs = params->advanfuncs;
	let options = params->options;
	let popmodel = &params->test.popmodel;
	var test = &params->test;
	let n = test->nparam;

	info(params->outstream, "deriv:\n");

	let f0 = params->best.result.objfn;

	var step2 = params->options->estimate.step_initial;
	var step1 = params->options->estimate.step_refine;
	forcount(k, n) {
		double x[OPENPMX_THETA_MAX + OPENPMX_OMEGA_MAX * OPENPMX_OMEGA_MAX + OPENPMX_SIGMA_MAX] = { };
		double xa[5] = { };
		double ya[5] = { };
		xa[2] = 0.;
		ya[2] = 0.;

		memset(x, 0, sizeof(x));

		x[k] = -1. * step2;
		encode_update(test, x);
		reset_eta(idata, params->besteta); 
		encode_evaluate(test, idata, advanfuncs, options);
		popmodel->result.neval = params->neval;
		params->neval += 1;
		popmodel->result.neval = params->neval;
		let fm2 = popmodel->result.objfn;
		print_model(params);
		xa[0] = x[k];
		ya[0] = fm2 - f0;
		x[k] = 0.;

		x[k] = -1. * step1;
		encode_update(test, x);
		reset_eta(idata, params->besteta); 
		encode_evaluate(test, idata, advanfuncs, options);
		params->neval += 1;
		popmodel->result.neval = params->neval;
		let fm1 = popmodel->result.objfn;
		print_model(params);
		xa[1] = x[k];
		ya[1] = fm1 - f0;
		x[k] = 0.;

		xa[2] = 0.;
		ya[2] = 0.;

		x[k] = 1. * step1;
		encode_update(test, x);
		reset_eta(idata, params->besteta); 
		encode_evaluate(test, idata, advanfuncs, options);
		params->neval += 1;
		popmodel->result.neval = params->neval;
		let fp1 = popmodel->result.objfn;
		print_model(params);
		xa[3] = x[k];
		ya[3] = fp1 - f0;
		x[k] = 0.;

		x[k] = 1. * step2;
		encode_update(test, x);
		reset_eta(idata, params->besteta); 
		encode_evaluate(test, idata, advanfuncs, options);
		params->neval += 1;
		popmodel->result.neval = params->neval;
		let fp2 = popmodel->result.objfn;
		print_model(params);
		xa[4] = x[k];
		ya[4] = fp2 - f0;
		x[k] = 0.;

		var interp = gsl_spline_alloc(gsl_interp_cspline, 5);
		var accelp = gsl_interp_accel_alloc();
		gsl_spline_init(interp, xa, ya, 5);

		let deriv = gsl_spline_eval_deriv(interp, 0., accelp);
		let deriv2 = gsl_spline_eval_deriv2(interp, 0., accelp);
		info(params->outstream, "param %i deriv %g %g (%g)\n", k, deriv, deriv2, deriv / deriv2);

		gsl_interp_accel_free(accelp);
		gsl_spline_free(interp);
	}
	/* TODO : add warnings here for decreases in objfn and zero gradient and second deriv */
	
	params->best.result.neval = params->neval;
}
#endif
	
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
	params->neval = 0; 
	if (!stabilize_model(params))
		warning(outstream, "initial evaluation not stable\n");

	/* make sure best has the objfn of what we just stabilized */
	params->best = params->test.popmodel;
	encode_offset(&params->test, &params->best); /* this probably does nothing */

	/* call the underlying advan to optimize and find the best imodel */
	let maxeval = options->estimate.maxeval;
	params->best.result.type = OBJFN_EVALUATE;
	if (maxeval > 1) {
		focei(params); /* TODO: this should be renamed I think */
		params->best.result.type = OBJFN_FINAL;
		params->best.result.neval = params->neval;

		/* at end we encode the best so far */
		encode_offset(&params->test, &params->best);
	}

	/* posthoc evaluation of derivates */
//	if (!options->estimate.posthoc.gradient.omit) 
//		evaluate_gradient(params);
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
	
	info(f2, "$OUTFILE\nOpenPMX %i.%i.%i %s\n", OPENPMX_VERSION_MAJOR, OPENPMX_VERSION_MINOR, OPENPMX_VERSION_RELEASE, OPENPMX_GITHASH);

#if defined(OPENPMX_PARALLEL_PTHREADS)
	char message[] = "pthread";
#elif defined(OPENPMX_PARALLEL_OPENMP)
	char message[] = "OpenMP";
#elif defined(OPENPMX_PARALLEL_SINGLETHREAD)
	char message[] = "single";
#else
#error no parallel processing advan defined
#endif
	info(f2, "config %s %i\n", message, options->nthread);

	info(f2, "data records %i used %i removed %i\n", advanfuncs->recordinfo.dataconfig->nrecords, idata->ndata, advanfuncs->recordinfo.dataconfig->nrecords - idata->ndata);
	info(f2, "data individuals %i observations %i\n", idata->nindivid, idata->nobs);
	info(f2, "data table offset %s\n", advanfuncs->recordinfo.dataconfig->_offset1 ? "true" : "false");

	advanfuncs->info(advanfuncs, stdout);
	if (f2)
		advanfuncs->info(advanfuncs, f2);

	if (outfile_type == OUTFILE_HEADER_RESAMPLE)
		info(f2, "resample seed %i\n", options->simulate.seed);
}

static STAGE2_PARAMS stage2_params(const char* filename,
								   IDATA* const idata,
								   const ADVANFUNCS* const advanfuncs,
								   POPMODEL* const popmodel,
								   const OPTIONS* const options)
{
	/* setup output logging */
	FILE* outstream = 0;
	FILE* extstream = 0;
	if (filename) {
		outstream = results_fopen(filename, OPENPMX_OUTFILE, "w");
		if (!outstream)
			fatal(0, "%s: could not open file \"%s%s\"\n", __func__, filename, OPENPMX_OUTFILE);
		extstream = results_fopen(filename, OPENPMX_EXTFILE, "w");
		if (!extstream)
			fatal(outstream, "%s: could not open file %s extension %s\n", __func__, filename, OPENPMX_EXTFILE);
	}

	let neta = idata->nindivid * idata->nomega;
	var params = (STAGE2_PARAMS) {
		.idata = idata,
		.advanfuncs = advanfuncs,
		.test = encode_init(popmodel),
		.options = options,
		.neval = 0,
		.best = *popmodel,
		.besteta = 0,
		.neta = neta,
		.begin = { }, 							/* set after initialization */
		.outstream = outstream,
		.extstream = extstream,
		.filename = filename,
	};
	clock_gettime(CLOCK_REALTIME, &params.begin);
	assert(idata->nindivid > 0);
	assert(idata->nomega > 0);

	params.best = params.test.popmodel; /* will set objfn to invalid */
	params.besteta = callocvar(double, neta);

	/* make sure everything is consistent */
	encode_offset(&params.test, &params.best);

	return params;
}

static void stage2_params_cleanup(STAGE2_PARAMS *params)
{
	if (params->outstream)
		fclose(params->outstream);
	if (params->extstream)
		fclose(params->extstream);
	free(params->besteta);
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
		message = focei(0);
		outfile_type = OUTFILE_HEADER_ESTIMATE;
	}

	/* setup the minimizer */
	var params = stage2_params(filename, idata, advanfuncs, popmodel, options);

	/* do some logging */
	/* start extfile and other headers, rest will be saved during iterations */
	outfile_header(params.outstream, advanfuncs, idata, options, outfile_type);
	var _offset1 = advanfuncs->recordinfo.dataconfig->_offset1;
	if (params.extstream) 
		extfile_header(params.extstream, &params.best, _offset1);
	if (maxeval > 1)
		info(params.outstream, "optim %s\n", message);

	/* do the checkout */
	idata_ineval(idata, true);
	idata_checkout(idata, advanfuncs, popmodel, options, params.outstream);

	/* actually do the stage 2 estimation */
	focei_popmodel_stage2(&params);
	*popmodel = params.best;

	/* update results to screen and log */
	let timestamp = get_timestamp(&params);
	popmodel_information(params.outstream, popmodel, timestamp);

	/* update results in tables */
	if (filename) {
		table_phi_idata(filename, idata, _offset1);
		if (params.extstream) {
			let ineval = idata_ineval(idata, false);
			extfile_trailer(params.extstream, &params.best, timestamp, ineval);
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

	var popmodel = popmodel_init(pmx);

	estimate_popmodel(pmx->filename,
					  &pstate->idata,
					  pstate->advanfuncs,
					  &popmodel,
					  &options);

	pmx_update_from_popmodel(pmx, &popmodel);
}

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

	var popmodel = popmodel_init(pmx);

	estimate_popmodel(pmx->filename,
					  &pstate->idata,
					  pstate->advanfuncs,
					  &popmodel,
					  &options);

	pmx_update_from_popmodel(pmx, &popmodel);
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

	var popmodel = popmodel_init(pmx);

	estimate_popmodel(pmx->filename,
					  &pstate->idata,
					  pstate->advanfuncs,
					  &popmodel,
					  &options);

	pmx_update_from_popmodel(pmx, &popmodel);
}


