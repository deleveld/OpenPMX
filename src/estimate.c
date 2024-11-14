/* 
 * This file is part of OpenPMX (https://github.com/deleveld/openpmx).
 * Copyright (c) 2022 Douglas Eleveld.
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
#include <values.h>

#include "openpmx.h"
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
		double term4 = omegainfo->omega_nonzero_lndet;
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

static double calculate_maxd(const int xlength,
							 const double x[static xlength])
{
	var d = 0.;
	if (x && xlength > 0) {
		forcount(j, xlength) {
			let delta = fabs(x[j]);
			if (d < delta)
				d = delta;
		}
	}
	return d;
}

typedef struct {
	IDATA* const idata;
	const ADVANFUNCS* const advanfuncs;
	ENCODE test;
	const OPTIONS* const options;

	int nfunc;
	POPMODEL best;
	double* besteta;

	struct timespec begin;
	FILE* outstream;
	const char* filename;
} STAGE2_PARAMS;

static void update_best_imodel(const int xlength,
							   const double x[static xlength],
							   const STAGE2_PARAMS* const params,
							   POPMODEL* const best)
{
	let idata = params->idata;
	let popmodel = &params->test.popmodel;

	const POPMODEL* improved_model = 0;
	let dobjfn = popmodel->result.objfn - best->result.objfn;
	if (dobjfn < 0.) {

		/* update the best estimation and its objective function and function evaluations so far */
		/* save the best eta so we can keep restarting there for speed and hopefully some consistancy */
		*best = *popmodel;
		improved_model = best;
		if (params->besteta) {
			let firstindivid = &idata->individ[0];
			memcpy(params->besteta, firstindivid->eta, idata->nindivid * idata->nomega * sizeof(double));
		}
	}

	/* update the user */
	let options = params->options;
	if (options->verbose && !options->brief)
		improved_model = popmodel;
	if (improved_model) {
		let maxd = calculate_maxd(xlength, x);

		struct timespec now;
		clock_gettime(CLOCK_REALTIME, &now);
		let runtime_s = timespec_time_difference(&params->begin, &now) / 1000.;
		char message[1024] = { 0 };
		sprintf(message, " dobjfn: %f", dobjfn);
		popmodel_eval_information(improved_model,
								  runtime_s,
								  params->filename,
								  options->verbose,
								  options->brief,
								  params->outstream,
								  xlength, x,
								  maxd, message);
	}
}

static void encode_evaluate(ENCODE* const test,
							IDATA* const idata,
							const ADVANFUNCS* const advanfuncs,
							const OPTIONS* const options)
							
{
	var popmodel = &test->popmodel;
	let omegainfo = &test->omegainfo;
	let nonzero = &omegainfo->nonzero;
	SCATTEROPTIONS scatteroptions = { 0 };
	scatteroptions.stage1_order = true;
	scatter_threads(idata, advanfuncs, popmodel, nonzero, options, &scatteroptions, stage1_thread);

	popmodel->result = (PMXRESULT) { .objfn = objfn(idata, omegainfo),
									 .type = OBJFN_CURRENT,
									 .nfunc = 0 };
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

	/* do the internal stage 1, start with eta from best run */
	assert(params->besteta);
	var firstindivid = &idata->individ[0];
	memcpy(firstindivid->eta, params->besteta, idata->nindivid * idata->nomega * sizeof(double));

	/* do the actual test, this sets the objfn */
	let advanfuncs = params->advanfuncs;
	let options = params->options;
	let popmodel = &params->test.popmodel;
	encode_evaluate(&params->test, idata, advanfuncs, options);
	params->nfunc += 1;
	popmodel->result.nfunc = params->nfunc;

	/* update best imodel and inform the user if we improve */
	update_best_imodel(_xlength, _x, params, &params->best);

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
	let step_size = params->options->estimate.optim.step_final;
	let n = params->test.nparam;
	var xx = gsl_vector_alloc(n);

	forcount(i, n) {
		let v = gsl_vector_get(x, i);
		volatile let h = step_size;
		volatile let above1 = v + h;
		volatile let above2 = v + 2. * h;
		volatile let below1 = v - h;
		volatile let below2 = v - 2. * h;

		gsl_vector_memcpy(xx, x);

		gsl_vector_set(xx, i, above2);
		volatile let f_above2 = foce_stage2_f(xx, _params);

		gsl_vector_set(xx, i, above1);
		volatile let f_above1 = foce_stage2_f(xx, _params);

		gsl_vector_set(xx, i, below1);
		volatile let f_below1 = foce_stage2_f(xx, _params);

		gsl_vector_set(xx, i, below2);
		volatile let f_below2 = foce_stage2_f(xx, _params);

//		let deriv = (f_above1 - f_below1) / (above1 - below1);	/* arbitrary method */

		/* http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/#noiserobust_2 */
		let deriv = (2. * (f_above1 - f_below1) + f_above2 - f_below2) / (8. * h);

		/* http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/lanczos-low-noise-differentiators/#lanczos_2 */
//		let deriv = ((f_above1 - f_below1) + 2. * (f_above2 - f_below2)) / (10. * h);

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
	let maxeval = options->estimate.optim.maxeval;
	let step_initial = options->estimate.optim.step_initial;
	let step_refine = options->estimate.optim.step_refine;
	let step_final = options->estimate.optim.step_final;

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
			if (params->nfunc >=  maxeval) 
				break;
		}
		gsl_multimin_fdfminimizer_free(s);

		let objfn = params->best.result.objfn;
		if (prevobjfn - objfn < 0.01)
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

	forcount(i, n)
		initial[i] = 0.;
	var neval = maxeval;
	var rhobeg = step_initial;
	var rhoend = step_refine;
	info(0, "initial: rho [%g, %g]\n", rhobeg, rhoend);
	var retcode = bobyqa(n, npt,
						 focei_stage2_evaluate_population_objfn, (void*)params,
						 initial, lower, upper,
						 rhobeg, rhoend,
						 iprint, neval, w);
	if (retcode != BOBYQA_SUCCESS)
		info(params->outstream, "initial BOBYQA error code %i\n", retcode);

	if (step_final < step_refine) {
		let best = &params->best;
		var lastobjfn = best->result.objfn;
		var dobjfn = DBL_MAX - lastobjfn;
		do {
			encode_offset(&params->test, &params->best);
			forcount(i, n)
				initial[i] = 0.;

			neval = maxeval - params->nfunc;
			rhobeg = step_refine;
			rhoend = step_final;
			info(params->outstream, "refine: rho [%g, %g]\n", rhobeg, rhoend);
			retcode = bobyqa(n, npt,
							 focei_stage2_evaluate_population_objfn, (void*)params,
							 initial, lower, upper,
							 rhobeg, rhoend,
							 iprint, neval, w);
			if (retcode != BOBYQA_SUCCESS)
				info(params->outstream, "refine BOBYQA error code %i\n", retcode);

			dobjfn = best->result.objfn - lastobjfn;
			info(params->outstream, "dobjfn: %f\n", dobjfn);
			lastobjfn = best->result.objfn;
		}
		while (dobjfn < -0.01);
	}

	free(w);
#endif

	free(lower);
	free(upper);
	free(initial);

	return 0;
}

static void print_model(STAGE2_PARAMS* params, const char* suffix)
{
	let popmodel = &params->test.popmodel;
	let options = params->options;
		
	struct timespec now;
	clock_gettime(CLOCK_REALTIME, &now);
	let runtime_s = timespec_time_difference(&params->begin, &now) / 1000.;
	popmodel_eval_information(popmodel,
							  runtime_s,
							  params->filename,
							  options->verbose,
							  options->brief,
							  params->outstream,
							  0, 0, 0,
							  suffix);
}

static bool stabilize_model(STAGE2_PARAMS* params)
{
	let idata = params->idata;
	let advanfuncs = params->advanfuncs;
	let options = params->options;
	
	var done = false;
	var niter = 0;
	while (!done) {
		let popmodel = &params->test.popmodel;
		let prev_objfn = popmodel->result.objfn;
		encode_evaluate(&params->test, idata, advanfuncs, options);
		params->nfunc += 1;
		++niter;
		popmodel->result.nfunc = params->nfunc;

		char message[1024] = { 0 };
		sprintf(message, " dobjfn: %f", popmodel->result.objfn - params->best.result.objfn);
		print_model(params, message);

		/* are we stable? */
		/* besteta we update later, individual eta values keep getting updated */
		if (fabs(popmodel->result.objfn - prev_objfn) < 0.01) 
			done = true;

		/* warn if we are not stable after 10 iterations */
		if (!done && niter >= 10)
			break;
	}
	return done;
}

static void evaluate_deriv(STAGE2_PARAMS* params)
{
	let idata = params->idata;
	let advanfuncs = params->advanfuncs;
	let options = params->options;

	let f0 = params->best.result.objfn;
	info(params->outstream, "deriv f0 %f\n", f0);

	var step = params->options->estimate.posthoc.deriv.step;
	var test = &params->test;
	forcount(k, test->nparam) {
		let popmodel = &params->test.popmodel;
		double x[OPENPMX_THETA_MAX + OPENPMX_OMEGA_MAX * OPENPMX_OMEGA_MAX + OPENPMX_SIGMA_MAX] = { 0 };

		x[k] = 1. * step;
		encode_update(test, x);
		encode_evaluate(test, idata, advanfuncs, options);
		params->nfunc += 1;
		let fp1 = popmodel->result.objfn;
		char message[1024] = { 0 };
		sprintf(message, " dobjfn: %f", fp1 - f0);
		print_model(params, message);
		x[k] = 0.;

		x[k] = -1. * step;
		encode_update(test, x);
		encode_evaluate(test, idata, advanfuncs, options);
		params->nfunc += 1;
		let fm1 = popmodel->result.objfn;
		sprintf(message, " dobjfn: %f", fm1 - f0);
		print_model(params, message);
		x[k] = 0.;

		x[k] = 2. * step;
		encode_update(test, x);
		encode_evaluate(test, idata, advanfuncs, options);
		params->nfunc += 1;
		let fp2 = popmodel->result.objfn;
		sprintf(message, " dobjfn: %f", fp2 - f0);
		print_model(params, message);
		x[k] = 0.;

		x[k] = -2. * step;
		encode_update(test, x);
		encode_evaluate(test, idata, advanfuncs, options);
		params->nfunc += 1;
		let fm2 = popmodel->result.objfn;
		sprintf(message, " dobjfn: %f", fm2 - f0);
		print_model(params, message);
		x[k] = 0.;

		x[k] = 3. * step;
		encode_update(test, x);
		encode_evaluate(test, idata, advanfuncs, options);
		params->nfunc += 1;
		let fp3 = popmodel->result.objfn;
		sprintf(message, " dobjfn: %f", fp3 - f0);
		print_model(params, message);
		x[k] = 0.;

		x[k] = -3. * step;
		encode_update(test, x);
		encode_evaluate(test, idata, advanfuncs, options);
		params->nfunc += 1;
		let fm3 = popmodel->result.objfn;
		sprintf(message, " dobjfn: %f", fm3 - f0);
		print_model(params, message);
		x[k] = 0.;

		/* http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/#noiserobust_2 */
		let deriv1 = (5.*(fp1 - fm1) + 4.*(fp2 - fm2) + fp3 - fm3) / (32. * step);

		/* https://en.wikipedia.org/wiki/Finite_difference_coefficient */
		let deriv2 = ((fp3 + fm3) + 2.*(fp2 - fm2) - (fp1 + fm1) - 4.*f0) / (16. * pow(step, 2));

		info(params->outstream, "param %i deriv1 %g deriv2 %g\n", k, deriv1, deriv2);
		
	}
	/* TODO : add warnings here for decreases in objfn and zero gradient and second deriv */
	
	params->best.result.nfunc = params->nfunc;
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
	if (!params->options->estimate.stage1.omit_icov_resample)
		idata_alloc_icovresample(idata);
	else
		idata_free_icovresample(idata);

	/* first run so we can set objective function and yhat. */
	/* this is the first evaluation, the encode_offset at the best model is
	 * already done. We may have to evaluate several times for a stable objfn */
	params->nfunc = 0; 
	if (!stabilize_model(params))
		warning(outstream, "initial evaluation not stable\n");

	/* after the iterations, we save the besteta which we use to start with later */
	memcpy(params->besteta, firstindivid->eta, idata->nindivid * idata->nomega * sizeof(double));
	/* we dont really need this since we dont change popmodel or best since starting
	 * best = *popmodel;
	 * encode_offset(&params->test, best); */

	/* call the underlying advan to optimize and find the best imodel */
	let maxeval = options->estimate.optim.maxeval;
	params->best.result.type = OBJFN_EVALUATE;
	if (maxeval > 1) {
		focei(params); /* TODO: this should be renamed I think */
		params->best.result.type = OBJFN_FINAL;
		params->best.result.nfunc = params->nfunc;

		/* at end ew encode the best so far */
		encode_offset(&params->test, &params->best);
		if (!stabilize_model(params))
			warning(outstream, "final evaluation not stable\n");
	}

	/* posthoc evaluation of derivates */
	if (!options->estimate.posthoc.deriv.omit) 
		evaluate_deriv(params);
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
	
	info(f2, "$OUTFILE\nOpenPMX %i.%i.%i\n", OPENPMX_VERSION_MAJOR, OPENPMX_VERSION_MINOR, OPENPMX_VERSION_RELEASE);

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
	info(f2, "config table offset: %s\n", options->_offset1 ? "true" : "false");

	info(f2, "data records: %i used %i total %i removed\n", idata->ndata, advanfuncs->recordinfo.dataconfig->nrecords, advanfuncs->recordinfo.dataconfig->nrecords - idata->ndata);
	info(f2, "data individuals: %i\n", idata->nindivid);
	info(f2, "data observations: %i\n", idata->nobs);

	advanfuncs->info(advanfuncs, stdout);
	if (f2)
		advanfuncs->info(advanfuncs, f2);

	if (outfile_type == OUTFILE_HEADER_RESAMPLE)
		info(f2, "resample seed: %i\n", options->simulate.seed);

/*
	if (outfile_type == OUTFILE_HEADER_EVALUATE ||
		outfile_type == OUTFILE_HEADER_ESTIMATE) {
		let stage1 = &options->estimate.stage1;
		info(f2, "stage1 gradient step: %g\n", stage1->gradient_step);
		info(f2, "stage1 step initial: %g\n", stage1->step_initial);
		info(f2, "stage1 step refine: %g\n", stage1->step_refine);
		info(f2, "stage1 stop final: %g\n", stage1->step_final);
		info(f2, "stage1 omit icov resample: %s\n", stage1->omit_icov_resample ? "true" : "false");
		info(f2, "stage1 max evaluations: %i\n", stage1->maxeval);
	}
	if (outfile_type == OUTFILE_HEADER_ESTIMATE) {
		let optim = &options->estimate.optim;
		info(f2, "optim step initial: %g\n", optim->step_initial);
		info(f2, "optim step refine: %g\n", optim->step_refine);
		info(f2, "optim step final: %g\n", optim->step_final);
		info(f2, "optim max evaluations: %i\n", optim->maxeval);
	}
*/
}

static STAGE2_PARAMS stage2_params(const char* filename,
								   IDATA* const idata,
								   const ADVANFUNCS* const advanfuncs,
								   POPMODEL* const popmodel,
								   const OPTIONS* const options)
{
	/* setup output logging */
	FILE* outstream = 0;
	if (filename) {
		outstream = results_fopen(filename, OPENPMX_OUTFILE, "w");
		if (!outstream)
			fatal(outstream, "%s: could not open file \"%s%s\"\n", __func__, filename, OPENPMX_OUTFILE);
	}

	var params = (STAGE2_PARAMS) {
		.idata = idata,
		.advanfuncs = advanfuncs,
		.test = encode_init(popmodel),
		.options = options,
		.nfunc = 0,
		.best = *popmodel,
		.besteta = 0,
		.begin = { 0 }, 							/* set after initialization */
		.outstream = outstream,
		.filename = filename,
	};
	clock_gettime(CLOCK_REALTIME, &params.begin);
	assert(idata->nindivid > 0);
	assert(idata->nomega > 0);

	params.best = params.test.popmodel; /* will set objfn to invalid */
	params.besteta = callocvar(double, idata->nindivid * idata->nomega);

	/* make sure everything is consistent */
	encode_offset(&params.test, &params.best);

	return params;
}

static void stage2_params_cleanup(STAGE2_PARAMS *params)
{
	if (params->outstream)
		fclose(params->outstream);
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
	let maxeval = options->estimate.optim.maxeval;
	popmodel->result.type = OBJFN_INITIAL;
	if (maxeval > 1) {
		message = focei(0);
		outfile_type = OUTFILE_HEADER_ESTIMATE;
	}

	/* setup the minimizer */
	var params = stage2_params(filename, idata, advanfuncs, popmodel, options);

	/* do some logging */
	/* start extfile and other headers, rest will be saved during iterations */
	outfile_header(params.outstream, advanfuncs, idata, options, outfile_type);
	var _offset1 = options->_offset1;
	if (filename)
		extfile_header(filename, &params.best, _offset1);
	if (maxeval > 1)
		info(params.outstream, "optim %s\n", message);

	/* do the checkout */
	idata_checkout(idata, advanfuncs, popmodel, options, params.outstream);

	/* actually do the stage 2 estimation */
	focei_popmodel_stage2(&params);

	/* update results */
	*popmodel = params.best;
	if (filename) {
		table_phi_idata(filename, idata, _offset1);
		extfile_trailer(filename, &params.best);
		if (!options->estimate.stage1.omit_icov_resample)
			table_icov_resample_idata(filename, idata, _offset1);
	}
	popmodel_information(params.outstream, popmodel); /* writes to stdout too */
	if (params.outstream) 
		popmodel_initcode(params.outstream, popmodel);

	/* cleanup */
	stage2_params_cleanup(&params);
}

void pmx_estimate(OPENPMX* pmx, ESTIMCONFIG* const estimate)
{
	pmxstate_ensure(pmx);
	var pstate = pmx->state;

	var options = options_init_from_pmx(pmx);
	if (estimate) {
		options.estimate = estimconfig_default(estimate);
		*estimate = options.estimate;
	}

	var popmodel = popmodel_init(pmx->theta, pmx->omega, pmx->sigma);
	popmodel.result.type = OBJFN_INITIAL;

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

	var options = options_init_from_pmx(pmx);
	if (stage1) {
		options.estimate.stage1 = stage1config_default(stage1);
		*stage1 = options.estimate.stage1;
	}
	options.estimate.optim.maxeval = 1;
	options.estimate.posthoc.deriv.omit = true;

	var popmodel = popmodel_init(pmx->theta, pmx->omega, pmx->sigma);
	popmodel.result.type = OBJFN_INITIAL;

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

	var options = options_init_from_pmx(pmx);
	if (estimate) {
		options.estimate = estimconfig_default(estimate);
		*estimate = options.estimate;
	}
	options.estimate.stage1.omit_icov_resample = true;
	options.estimate.optim.step_initial = 1.;
	options.estimate.optim.step_refine = 0.1;
	options.estimate.optim.step_final = 0.01;
	options.estimate.posthoc.deriv.omit = true;

	var popmodel = popmodel_init(pmx->theta, pmx->omega, pmx->sigma);
	popmodel.result.type = OBJFN_INITIAL;

	estimate_popmodel(pmx->filename,
					  &pstate->idata,
					  pstate->advanfuncs,
					  &popmodel,
					  &options);

	pmx_update_from_popmodel(pmx, &popmodel);
}


