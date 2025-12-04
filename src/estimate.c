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
#include "options.h"
#include "pmxstate.h"
#include "advan/advan.h"
#include "utils/c22.h"
#include "utils/various.h"

#include "openpmx_compile_options.h"

/*--------------------------------------------------------------------*/
/* different optimizers */

#define OPTIMIZER_OUTER_BOBYQA
//#define OPTIMIZER_OUTER_GSL_NELDERMEAD
//#define OPTIMIZER_OUTER_GSL_BFGS
// TODO: Not working, sigsev when outer is libprima as well!
// The data pointer getting passed to inner is invalid
//#define OPTIMIZER_OUTER_LIBPRIMA

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
		let ineval = idata_ineval(idata, false);
		var outstream = (options->estimate.progress) ? (params->outstream) : 0;
		var extstream = (options->estimate.progress) ? (params->extstream) : 0;
		let runtime_s = get_timestamp(params);
		popmodel_eval_information(improved_model,
								  runtime_s,
								  ineval,
								  options->estimate.details || options->estimate.verbose,
								  outstream,
								  extstream,
								  0);
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
	SCATTEROPTIONS scatteroptions = { };
	scatteroptions.stage1_order = true;
	scatter_threads(idata, advanfuncs, popmodel, nonzero, options, &scatteroptions, stage1_thread);

	popmodel->result = (PMXRESULT) {
		.objfn = idata_objfn(idata, omegainfo->omega_nonzero_lndet),
		.type = OBJFN_CURRENT,
		.nparam = 0,
		.neval = 0
	};
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
	idata_reset_eta(idata, params->besteta);
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

#ifdef OPTIMIZER_OUTER_LIBPRIMA
//#include "prima.h"
static void outer_fun(const double x[], double* const f, const void* data)
{
//	info(0, "Outer data %p\n", data);
//	fflush(stdout);

	let params = (STAGE2_PARAMS*) data;
	let objfn = focei_stage2_evaluate_population_objfn(params->test.nparam,
													   x, 
													   data);
    *f = objfn;
}
#endif

static const char* focei(STAGE2_PARAMS* const params)
{
	if (!params)
		return "FOCEI BOBYQA";

	let options = params->options;

	let n = params->test.nparam;
	var initial = mallocvar(double, n);
	var lower = mallocvar(double, n);
	var upper = mallocvar(double, n);
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

#ifdef OPTIMIZER_OUTER_LIBPRIMA
	let best = &params->best;
	var lastobjfn = best->result.objfn;
	forcount(i, n)
		initial[i] = 0.;
	var neval = maxeval;
	var rhobeg = step_initial;
	var rhoend = step_refine;
	info(params->outstream, "optim rho %g %g\n", rhobeg, rhoend);
	
	prima_problem_t problem;
	prima_init_problem(&problem, n);
	problem.x0 = initial;
	problem.xl = lower;
	problem.xu = upper;
	problem.calfun = outer_fun;

	// Set up the options
	prima_options_t poptions;
	prima_init_options(&poptions);
	poptions.rhobeg = rhobeg;
	poptions.rhoend = rhoend;
	poptions.maxfun = neval;
	poptions.data = (void*)params;
	poptions.callback = 0;

	// initial rough estimation
#define PRIMA_METHOD PRIMA_BOBYQA
#define PRIMA_METHOD PRIMA_LINCOA
	prima_result_t result;
	prima_rc_t retcode = prima_minimize(PRIMA_METHOD, problem, poptions, &result);
    prima_free_result(&result);

//	if (retcode != PRIMA_RESULT_INITIALIZED)
//		warning(0, "initial BOBYQA(libprima) error code %i\n", retcode);

	var runtime_s = get_timestamp(params);
	info(params->outstream, "time %.3f neval %i objfn %f\n", runtime_s, params->neval, best->result.objfn);
	lastobjfn = best->result.objfn;

	// refine estimation
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

			poptions.rhobeg = rhobeg;
			poptions.rhoend = rhoend;
			poptions.maxfun = neval;
			retcode = prima_minimize(PRIMA_METHOD, problem, poptions, &result);
			prima_free_result(&result);

//			if (retcode != PRIMA_RESULT_INITIALIZED)
//				warning(0, "refine BOBYQA(libprima) error code %i\n", retcode);

			dobjfn = best->result.objfn - lastobjfn;
			var runtime_s = get_timestamp(params);
			info(params->outstream, "time %.3f neval %i objfn %f\n", runtime_s, params->neval, best->result.objfn);
			lastobjfn = best->result.objfn;
		}
		while (dobjfn < -1.*fabs(options->estimate.dobjfn));
	}
	info(params->outstream, "optim rho %g\n", rhoend);

    prima_free_result(&result);
#endif

	free(lower);
	free(upper);
	free(initial);

	return 0;
}

static void estimate_print_model(STAGE2_PARAMS* params)
{
	let idata = params->idata;
	let popmodel = &params->test.popmodel;
	let options = params->options;

//	char message[128] = "";
//	if (popmodel->result.objfn != DBL_MAX && params->best.result.objfn != DBL_MAX)
//		sprintf(message, " objfn %f", popmodel->result.objfn);

	let ineval = idata_ineval(idata, false);
	var outstream = (options->estimate.progress) ? (params->outstream) : 0;
	var extstream = (options->estimate.progress) ? (params->extstream) : 0;
	let runtime_s = get_timestamp(params);
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
	idata_reset_eta(idata, params->besteta);

	/* very first evaluation */
	encode_evaluate(&params->test, idata, advanfuncs, options);
	params->neval += 1;
	var timestamp = get_timestamp(params);
	info(params->outstream, "time %.3f neval %i objfn %f\n", timestamp, params->neval, popmodel->result.objfn);

	var done = false;
	var niter = 0;
	let maxiter = 10;
	while (!done) {
		let lastobjfn = popmodel->result.objfn;
		encode_evaluate(&params->test, idata, advanfuncs, options);
		params->neval += 1;
		++niter;
		popmodel->result.neval = params->neval;

		estimate_print_model(params);

		/* are we stable? */
		/* besteta we update later, individual eta values keep getting updated */
		/* warn if we are not stable after 10 iterations */
		let dobjfn = popmodel->result.objfn - lastobjfn;
		if (fabs(dobjfn) < 0.01) 
			done = true;
		else if (niter >= maxiter) 
			break;
		params->best = *popmodel;
	}
	info(params->outstream, "stabilize iter %i\n", niter);
	if (niter == maxiter) 
		warning(params->outstream, "stabilize failed\n");
	timestamp = get_timestamp(params);
	info(params->outstream, "time %.3f neval %i objfn %f\n", timestamp, params->neval, popmodel->result.objfn);

	var firstindivid = &idata->individ[0];
	memcpy(params->besteta, firstindivid->eta, idata->nindivid * idata->nomega * sizeof(double));

	return done;
}

#if 0
#include <gsl/gsl_spline.h>

typedef struct {
	double* point;
	int dimnum;
	double objfn;
} COVPOINTS;

#include "linalg.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

void print_matrix(const gsl_matrix* const p)
{
	if (!p) {
		printf("(null)\n");
		return;
	}
	printf("matrix(%i,%i)\n", (int)p->size1, (int)p->size2);
	forcount(i, p->size1) {
		forcount(j, p->size2)
			printf("%13e ", gsl_matrix_get(p, i, j));
		printf("\n");
	}
}

static void reconstruct(const double* _x, const int _xlength, const int n, const int nchol, gsl_matrix* chol2, double* mean2, double* offset)
{
	var i = 0;
	*offset = _x[i++];
	forcount(ii, n)
		mean2[ii] = _x[i++];
	gsl_matrix_set_zero(chol2);
	var r = 0;
	var c = 0;
	assert(nchol == n);
	forcount(ii, nchol) {
		var v = _x[i++];
		r = ii;
		c = ii;
		if (r == c)
			v = exp(v);
		gsl_matrix_set(chol2, r, c, v);
	}
	assert(i == _xlength);
} 

typedef struct {
	const COVPOINTS* const covpoints;
	const int ncovpoints;
	const int n;
	const int nchol;
	double bestf;
} COVINFO;

static int foce_gradient_f(const gsl_vector* const x, void *data, gsl_vector * f)
{
	let covinfo = (COVINFO*) data;
	let n = covinfo->n;
	let nchol = covinfo->nchol;

	var off2 = 0.;
	var mean2 = callocvar(double, n);
	var chol2 = gsl_matrix_alloc(n, n);
	reconstruct(x->data, x->size, n, nchol, chol2, mean2, &off2);

	let covpoints = covinfo->covpoints;
	var ssqerr = 0.;
	forcount(i, covinfo->ncovpoints) {
		let p = covpoints[i].point;
		double adj[n];
		forcount(k, n)
			adj[k] = p[k] - mean2[k];
		let pobjfn = sample_min2ll_from_cholesky(adj, chol2) + off2;
		let err = covpoints[i].objfn - pobjfn;
		
		gsl_vector_set(f, i, err);
		ssqerr += pow(err, 2);
	}
	info(0, "ssqerr=%5f\n", ssqerr);
	
	gsl_matrix_free(chol2);
	free(mean2);
	
	return GSL_SUCCESS;
}

static int foce_gradient_df(const gsl_vector * x, void *data, gsl_matrix * J)
{
	let n = J->size1;
	let p = x->size;
	gsl_vector *xx = gsl_vector_alloc(p);
	gsl_vector *f_plus_h = gsl_vector_alloc(n);
	gsl_vector *f_minus_h = gsl_vector_alloc(n);

	forcount(j, p) {
		let step_size = 0.1;

		let v = gsl_vector_get(x, j);
		let h = step_size;
		let above = v + h;
		let below = v - h;
		let two_times_h = above - below;

		gsl_vector_memcpy(xx, x);

		gsl_vector_set(xx, j, above);
		foce_gradient_f(xx, data, f_plus_h);

		gsl_vector_set(xx, j, below);
		foce_gradient_f(xx, data, f_minus_h);

		forcount(i, n) {
			let upper = gsl_vector_get(f_plus_h, i);
			let lower = gsl_vector_get(f_minus_h, i);
			let deriv = (upper - lower) / (two_times_h);
			gsl_matrix_set(J, i, j, deriv);
		}
	}
	gsl_vector_free(xx);
	gsl_vector_free(f_plus_h);
	gsl_vector_free(f_minus_h);
	return GSL_SUCCESS;
}

#include <gsl/gsl_multifit_nlinear.h>
static void foce_gradient_callback(const size_t iter, void * params, const gsl_multifit_nlinear_workspace * w)
{
	(void)w;
	int* niter = (int*) params;
	*niter = iter;
}

static void minimize_covpoints(const COVPOINTS* const covpoints, const int ncovpoints, const int n, const int nchol, double* x, double* mean2, gsl_matrix* chol2, double* off2)
{
	var covinfo = (COVINFO) {
		.covpoints = covpoints,
		.ncovpoints = ncovpoints,
		.n = n,
		.nchol = nchol,
		.bestf = DBL_MAX,
	};
	let ndim = 1 + n + nchol;

	var fdf = (gsl_multifit_nlinear_fdf) {
		.f = foce_gradient_f,
//		.df = foce_gradient_df, /* TODO: Can this be efficently done with step-size equal to stddev of observation? */
		.df = NULL, /* TODO: Can this be efficently done with step-size equal to stddev of observation? */
		.fvv = NULL,
		.n = (size_t)ncovpoints,
		.p = (size_t)ndim,
		.params = (void*)&covinfo,
		.nevalf = 0,
		.nevaldf = 0,
		.nevalfvv = 0,
	};
	
	var fdf_params = gsl_multifit_nlinear_default_parameters();
//	fdf_params.trs =  gsl_multifit_nlinear_trs_lm;
	fdf_params.trs =  gsl_multifit_nlinear_trs_lmaccel;
//	fdf_params.trs = gsl_multifit_nlinear_trs_dogleg;
//	fdf_params.trs = gsl_multifit_nlinear_trs_ddogleg;
//	fdf_params.trs = gsl_multifit_nlinear_trs_subspace2D;

	fdf_params.scale = gsl_multifit_nlinear_scale_more;
//	fdf_params.scale = gsl_multifit_nlinear_scale_levenberg;
//	fdf_params.scale = gsl_multifit_nlinear_scale_marquardt;

	fdf_params.solver = gsl_multifit_nlinear_solver_qr;
//	fdf_params.solver = gsl_multifit_nlinear_solver_cholesky;
//	fdf_params.solver = gsl_multifit_nlinear_solver_svd;

//	fdf_params.fdtype = GSL_MULTIFIT_NLINEAR_FWDIFF;
//	fdf_params.fdtype = GSL_MULTIFIT_NLINEAR_CTRDIFF;

	fdf_params.h_df = 0.1;
//		fdf_params.h_fvv = 0.1; // estimoptions->stage1_step_size;

	fdf_params.avmax = 0.4;

	let xtol = 1e-3;
	let gtol = pow(GSL_DBL_EPSILON, 1./3.); 	/* based on https://www.gnu.org/software/gsl/manual/html_node/Nonlinear-Least_002dSquares-Testing-for-Convergence.html#Nonlinear-Least_002dSquares-Testing-for-Convergence */
	let ftol = 0.;
	let max_iter = 1000;

	/* initialize solver with starting point and weights */
	var niter = 0;
	var w = gsl_multifit_nlinear_alloc(gsl_multifit_nlinear_trust, &fdf_params, ncovpoints, ndim);
	let initial = gsl_vector_const_view_array(x, ndim);
	gsl_multifit_nlinear_init(&initial.vector, &fdf, w);
	var info = -1;
	gsl_multifit_nlinear_driver(max_iter, 
		xtol, gtol, ftol,
		foce_gradient_callback, &niter,
		&info, w);

	/* final results */
	let xfinal = gsl_multifit_nlinear_position(w);
	forcount(i, ndim) 
		x[i] = gsl_vector_get(xfinal, i);
	reconstruct(x, ndim, n, nchol, chol2, mean2, off2);

	gsl_multifit_nlinear_free(w);
}

#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>

static void evaluate_gradient(STAGE2_PARAMS* params)
{
	return;
	
	/* http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/#noiserobust_2 */
	
	let idata = params->idata;
	let advanfuncs = params->advanfuncs;
	let options = params->options;
	let popmodel = &params->test.popmodel;
	var test = &params->test;
	let n = test->nparam;
	var outstream = params->outstream;
	var extstream = params->extstream;

	info(outstream, "deriv:\n");

	let baseobjfn = params->best.result.objfn;

	VECTOR(COVPOINTS) covpoints = { };
	var basepoint = (COVPOINTS) {
		.point = callocvar(double, n),
		.dimnum = 999,
		.objfn = baseobjfn
	};
	vector_append(covpoints, basepoint);
	
	let step0 = 0.02;
	double stepsize[n];
	forcount(i, n)
		stepsize[i] = 0.;
	forcount(k, n) {
		info(outstream, "paramater %i -------------------------\n", k);
		VECTOR(COVPOINTS) dirpoints = { };

		/* lower bound both sized of delta objfn */
		var step = -1. * step0;
		var niter = 0;
		var delta = 0.;
		while (delta < 10. && niter++ < 15) {
			double x[OPENPMX_THETA_MAX + OPENPMX_OMEGA_MAX * OPENPMX_OMEGA_MAX + OPENPMX_SIGMA_MAX] = { };
			x[k] = step;
			encode_update(test, x);
			idata_reset_eta(idata, params->besteta); 
			encode_evaluate(test, idata, advanfuncs, options);
			estimate_print_model(params);
			params->neval += 1;
			popmodel->result.neval = params->neval;
			delta = popmodel->result.objfn - baseobjfn;
			var thispoint = (COVPOINTS) {
				.point = callocvar(double, n),
				.dimnum = k,
				.objfn = popmodel->result.objfn,
			};
			thispoint.point[k] = x[k];
			vector_append(dirpoints, thispoint); 

			step *= 2.;
		}
		stepsize[k] = step / 2;

		/* upper bound both sized of delta objfn */
		step = 1. * step0;
		niter = 0;
		delta = 0.;
		while (delta < 10. && niter++ < 15) {
			double x[OPENPMX_THETA_MAX + OPENPMX_OMEGA_MAX * OPENPMX_OMEGA_MAX + OPENPMX_SIGMA_MAX] = { };

			x[k] = step;
			encode_update(test, x);
			idata_reset_eta(idata, params->besteta); 
			encode_evaluate(test, idata, advanfuncs, options);
			estimate_print_model(params);
			params->neval += 1;
			popmodel->result.neval = params->neval;
			delta = popmodel->result.objfn - baseobjfn;
			var thispoint = (COVPOINTS) {
				.point = callocvar(double, n),
				.dimnum = k,
				.objfn = popmodel->result.objfn,
			};
			thispoint.point[k] = x[k];
			vector_append(dirpoints, thispoint); 

			step *= 2.;
		}
		stepsize[k] = (stepsize[k] + step / 2.) / 2.;

		let npoints = vector_size(dirpoints);
		vector_appendn(covpoints, dirpoints.rawptr, npoints);
		vector_free(dirpoints);
	}

	var m = gsl_matrix_alloc(n, n);
	gsl_matrix_set_zero(m);
	forcount(j, n)
		gsl_matrix_set(m, j, j, 0.1);
	var chol = gsl_matrix_alloc(n, n);
	gsl_matrix_memcpy(chol, m);
	gsl_linalg_cholesky_decomp1(chol);

	print_matrix(m);
	print_matrix(chol);

//	let nchol = n * (n + 1) / 2;
	let nchol = n;
	let ndim = 1 + n + nchol;
	var offset = baseobjfn;
	var x = callocvar(double, ndim);
	var r = 0;
	var c = 0;
	var i = 0;
	x[i++] = offset;
	forcount(ii, n)			/* mean */
		x[i++] = 0.;
	forcount(ii, nchol) {	/* cholesky */
		r = ii;
		c = ii;
		var v = gsl_matrix_get(chol, r, c);
		if (r == c)
			v = log(v);
		x[i++] = v;
	}
	assert(i == 1 + n + nchol);
	assert(i == ndim);

	var off2 = 0.;
	var mean2 = callocvar(double, n);
	var chol2 = gsl_matrix_alloc(n, n);
	reconstruct(x, ndim, n, nchol, chol2, mean2, &off2);
	
	var f0 = 0.;
	forvector(i, covpoints) {
		let p = covpoints.ptr[i].point;
		double adj[n];
		forcount(k, n)
			adj[k] = p[k] - mean2[k];
		let pobjfn = sample_min2ll_from_cholesky(adj, chol2) + off2;
		let err = covpoints.ptr[i].objfn - pobjfn;
		f0 += (err * err);
		info(outstream, "dim %i objfn=%g pobjfn=%g err=%f\n",
			covpoints.ptr[i].dimnum,
			covpoints.ptr[i].objfn,
			pobjfn, err);  
	}
	info(0, "initial f0 = %g\n", f0);
	
	minimize_covpoints(covpoints.ptr,
		vector_size(covpoints),
		n, nchol,
		x, mean2, chol2, &off2);

	let npoints = 100;
	forcount(i, npoints) {
		static unsigned long int seed = 0;
		static gsl_rng* rng;
		if (seed == 0) {
			rng = gsl_rng_alloc(gsl_rng_mt19937);
			seed = 200501041406;
			gsl_rng_set(rng, seed);
		}
		var v = gsl_vector_alloc(n);
		forcount(j, n) {
			var x = gsl_ran_gaussian(rng, 0.5);
			if (i == npoints - 1)
				x = 0.;
			gsl_vector_set(v, j, x);
		}
		var s = gsl_vector_alloc(n);
		gsl_blas_dgemv(CblasNoTrans, 1., chol2, v, 0., s);
		var newsample = mallocvar(double, n);
		forcount(j, n) {
			let x = gsl_vector_get(s, j);
			newsample[j] = x + mean2[j]; 
		}
		gsl_vector_free(s);
		gsl_vector_free(v);

		encode_update(test, newsample);
		idata_reset_eta(idata, params->besteta); 
		encode_evaluate(test, idata, advanfuncs, options);
		estimate_print_model(params);
		params->neval += 1;
		popmodel->result.neval = params->neval;
		var newpoint = (COVPOINTS) {
			.point = callocvar(double, n),
			.dimnum = -1,
			.objfn = popmodel->result.objfn,
		};
		forcount(k, n)
			newpoint.point[k] = newsample[k];
		free(newsample);
		if (i == npoints - 1)
			newpoint.dimnum = -2;
		vector_append(covpoints, newpoint);

		info(0, "... solved this point\n");
		let timestamp = get_timestamp(params);
		popmodel_information(outstream, popmodel, timestamp);

		info(0, "... remove bad points\n");
		var minobjfn = DBL_MAX;
		forvector(j, covpoints) {
			var v = &covpoints.ptr[j];
			if (v->objfn < minobjfn)
				minobjfn = v->objfn;
		}
		printf("minobjfn = %f\n", minobjfn);
		forvector(j, covpoints) {
			var v = &covpoints.ptr[j];
			if (vector_size(covpoints) <= n + nchol + 20)
				break;
			if (v->objfn - minobjfn > 10.) {
				printf("remove = %i\n", j);
				vector_remove(covpoints, j, 1);
				--j;
			}
		}

		info(0, "... minimize\n");
		minimize_covpoints(covpoints.ptr,
			vector_size(covpoints),
			n, nchol,
			x, mean2, chol2, &off2);
			
		info(0, "... points during minimize\n");
		forvector(j, covpoints) {
			let p = covpoints.ptr[j].point;
			double adj[n];
			forcount(k, n)
				adj[k] = p[k] - mean2[k];
			let pobjfn = sample_min2ll_from_cholesky(adj, chol2) + off2;
			let err = covpoints.ptr[j].objfn - pobjfn;
			info(outstream, "%i dim %i objfn=%g pobjfn=%g err=%f\n",
				j,
				covpoints.ptr[j].dimnum,
				covpoints.ptr[j].objfn,
				pobjfn, 
				err);  
		}
	}
	info(0, "... points after minimize\n");
	forvector(j, covpoints) {
		let p = covpoints.ptr[j].point;
		double adj[n];
		forcount(k, n)
			adj[k] = p[k] - mean2[k];
		let pobjfn = sample_min2ll_from_cholesky(adj, chol2) + off2;
		let err = covpoints.ptr[j].objfn - pobjfn;
		info(outstream, "%i dim %i objfn=%g pobjfn=%g err=%f\n",
			j,
			covpoints.ptr[j].dimnum,
			covpoints.ptr[j].objfn,
			pobjfn, 
			err);  
	}
	
	info(0, "offset %f\nmean", off2);
	forcount(i, n)
		info(0, " %g", mean2[i]);
	info(0, "\ncov\n");
	gsl_blas_dgemm(CblasNoTrans, CblasTrans,
				   1.0, chol2, chol2,
				   0.0, m);
	print_matrix(m);
	
	var s = fopen("o.dat", "w");
	assert(s);
	fprintf(s, "DIMN\tOBJFN\tPOBJFN");
	forcount(j, n)
		fprintf(s, "\tX%i", j);
	fprintf(s, "\n");
	forvector(i, covpoints) {
		let p = covpoints.ptr[i].point;
		double adj[n];
		forcount(k, n)
			adj[k] = p[k] - mean2[k];
		let pobjfn = sample_min2ll_from_cholesky(adj, chol2) + off2;
		let err = covpoints.ptr[i].objfn - pobjfn;
		fprintf(s, "%i\t%f\t%f", 
			covpoints.ptr[i].dimnum,
			covpoints.ptr[i].objfn, 
			pobjfn);  
		forcount(j, n)
			fprintf(s, "\t%f", p[j]);
		fprintf(s, "\n");
	}
	fclose(s);
	
	s = fopen("o.R", "w");
	assert(s);
	fprintf(s, 
		"postscript(file=\"o.ps\", horizontal=TRUE)\n"
		"f <- read.table(file=\"o.dat\", header=TRUE)\n"
		"print(head(f))\n"
		"print(tail(f))\n"
		"plot(x=f$OBJFN,y=f$POBJFN)\n"
		"fsel<-f$DIMN==-1\n"
		"points(x=f$OBJFN[fsel],y=f$POBJFN[fsel], pch=19, cex=2, col=\"lightblue\")\n"
		"fsel<-f$DIMN==-2\n"
		"points(x=f$OBJFN[fsel],y=f$POBJFN[fsel], pch=19, cex=2, col=\"lightgreen\")\n"
		"par(mfrow=c(2,3))\n"
		"\n");
	forcount(i, n) {
		fprintf(s, 
			"fmin<-min(f$OBJFN)\n"
			"fsel<-f$DIMN==%i|f$DIMN==999\n"
			"ff<-f[fsel,]\n"
			"ylim<-range(fmin, fmin+10)\n"
			"plot(x=ff$X%i,y=ff$OBJFN,ylim=ylim)\n"
			"title(\"Dim %i\")\n"
			"points(x=ff$X%i,y=ff$OBJFN, pch=2)\n"
			"xextra<-f[f$DIMN==-1,]\n"
			"points(x=xextra$X%i,y=xextra$OBJFN, pch=19, cex=2, col=\"lightblue\")\n"
			"xextra<-f[f$DIMN==-2,]\n"
			"points(x=xextra$X%i,y=xextra$OBJFN, pch=19, cex=2, col=\"lightgreen\")\n",
			i, i, i, i, i, i);
	} 
	fclose(s);
	system("Rscript o.R");
	exit(8);

	encode_update(test, mean2);
	idata_reset_eta(idata, params->besteta); 
	encode_evaluate(test, idata, advanfuncs, options);
	estimate_print_model(params);
	params->neval += 1;
	popmodel->result.neval = params->neval;
	let timestamp = get_timestamp(params);
	popmodel_information(outstream, popmodel, timestamp);
	

	free(x);
	free(mean2);
	gsl_matrix_free(chol2);
	gsl_matrix_free(chol);
	gsl_matrix_free(m);
	
	forvector(i, covpoints)
		free(covpoints.rawptr[i].point);
	vector_free(covpoints);
	
	params->best.result.neval = params->neval;
	/* TODO : add warnings here for decreases in objfn and zero gradient and second deriv */
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
		params->best.result.nparam = params->test.nparam;
		params->best.result.neval = params->neval;

		/* at end we encode the best so far */
		encode_offset(&params->test, &params->best);
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
	popmodel->result.nparam = params.test.nparam;

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
	
	/* posthoc evaluation of derivates */
#if 0
//	if (!options->estimate.posthoc.gradient.omit) 
		evaluate_gradient(&params);
#endif
	
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


