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

#include <string.h>
#include <assert.h>
#include <float.h>
#include <math.h>

#include "stage1.h"
#include "ievaluate.h"
#include "nonzero.h"
#include "linalg.h"
#include "print.h"
#include "advan/advan.h"
#include "utils/c22.h"
#include "utils/various.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

/*--------------------------------------------------------------------*/
/* different optimizers */
#define OPTIMIZER_INNER_BOBYQA

/*--------------------------------------------------------------------*/

typedef struct {
	double* const testeta;
	const NONZERO* const nonzero;
	int* const ineval;
	const STAGE1CONFIG* const stage1;
	const IEVALUATE_ARGS ievaluate_args;
	double* const eval_msec;
	const int nobs;
} STAGE1_PARAMS;

static double stage1_evaluate_individual_iobjfn(const long int nreta,
												const double* reta,
												void* const data)
{
	let stage1_params = (const STAGE1_PARAMS * const) data;
	let ievaluate_args = &stage1_params->ievaluate_args;

	/* reta and nreta are the variable in the reduced eta/omega matricies */
	/* have to un-reduce the eta values because user code needs the full eta */
	/* the user code can see the testeta because ievaluate_args->popparam->eta
	 * points to the same place */
	assert(nreta == stage1_params->nonzero->n);
	unreduce_eta(stage1_params->testeta, reta, stage1_params->nonzero);

	/* functions for Bae and Yim objective function Term 1 and Term 2 */
	struct timespec t3;
	clock_gettime(CLOCK_MONOTONIC, &t3);
	let objfn_term1_term2 = individual_fasteval(ievaluate_args);
	timespec_duration(&t3, stage1_params->eval_msec);
	*(stage1_params->ineval) += 1;

	/* functions for Bae and Yim objective function Term 3 */
	var term3 = sample_min2ll(nreta, reta, stage1_params->nonzero);

	let iobjfn = objfn_term1_term2 + term3;
	return iobjfn;
}

#ifdef OPTIMIZER_INNER_BOBYQA
#include "bobyqa/bobyqa.h"
#endif

/* NOTE: this function must be thread safe and only touch individual data */
static bool estimate_individual_posthoc_eta(double reta[static OPENPMX_OMEGA_MAX],
											const STAGE1_PARAMS* const stage1_params)
{
	var nonzero = stage1_params->nonzero;
	let ievaluate_args = &stage1_params->ievaluate_args;
	let popparam = &ievaluate_args->popparam;
	let stage1 = stage1_params->stage1;

	/* detect special first case which could probably need more iterations */
	bool all_eta_zero = true;

	/* extract previous reduced eta as start point */
	assert(stage1_params->testeta == popparam->eta);
	reduce_eta(reta, popparam->eta, nonzero);
	double lower[OPENPMX_OMEGA_MAX];
	double upper[OPENPMX_OMEGA_MAX];
	let nreta = nonzero->n;
	forcount(i, nreta) {
		lower[i] = -1e6; // -DBL_MAX;
		upper[i] = 1e6; // DBL_MAX;
		if (reta[i] != 0.)
			all_eta_zero = false;
	}

	/* extra initial optimization when starting from 0 */
	let n = nreta;
	assert(n < OPENPMX_OMEGA_MAX);
	let neval = stage1->maxeval;
	
#ifdef OPTIMIZER_INNER_BOBYQA
	var rhobeg = stage1->step_initial;
	var rhoend = stage1->step_refine;
	int retcode;
	let iprint = 0;
	let npt = 2*n+1;				/* recommended */
//	let npt_min = n+2;				/* minimum */
//	let npt_max = (n+1)*(n+2)/2;	/* maximum */

	let wsize = (npt+5)*(npt+n)+3*n*(n+5)/2 + 10; 	/* a little bit extra room to be sure */
	var w = mallocvar(double, wsize);
	if (all_eta_zero) {
		retcode = bobyqa(n, npt,
						 stage1_evaluate_individual_iobjfn,
						 (void*)stage1_params,
						 reta, lower, upper,
						 rhobeg, rhoend,
						 iprint, neval, w);
		if (retcode != BOBYQA_SUCCESS) {
			let recordinfo = &stage1_params->ievaluate_args.advanfuncs->recordinfo;
			let id = RECORDINFO_ID(recordinfo, stage1_params->ievaluate_args.record);
			warning(0, "BOBYQA error %i: ID %f eta (initial) not successful\n", retcode, id);
		}
	}

	/* regular refinement optimization */
	/* we dont need to realloc because the previous memory in w is enough */
	rhobeg = stage1->step_refine;
	rhoend = stage1->step_final;
	retcode = bobyqa(n, npt,
					 stage1_evaluate_individual_iobjfn, (void*)stage1_params,
					 reta, lower, upper,
					 rhobeg, rhoend,
					 iprint, neval, w);
	if (retcode != BOBYQA_SUCCESS) {
		let recordinfo = &stage1_params->ievaluate_args.advanfuncs->recordinfo;
		let id = RECORDINFO_ID(recordinfo, stage1_params->ievaluate_args.record);
		warning(0, "BOBYQA error %i: ID %f eta (refine) not successful\n", retcode, id);
	}
	free(w);
#endif

	/* final results, i.e. the etas of the reduced matrix, get expanded */
	assert(stage1_params->testeta == popparam->eta);
	unreduce_eta(stage1_params->testeta, reta, nonzero);
	
	return all_eta_zero;
}

static void write_icov_from_reduced(double* icov, 
									const NONZERO* const nonzero,
									const gsl_matrix* const reducedcov,
									const int nomega)
{
	let nreta = reducedcov->size1;
	let rowcol = nonzero->rowcol;
	forcount(i, nreta) {
		forcount(j, nreta) {
			let row = rowcol[i];
			let col = rowcol[j];
			let v = gsl_matrix_get(reducedcov, i, j);
			icov[row * nomega + col] = v;
		}
	}
}

/* calculate individual variance covariance matrix
 * yhatvar must be non-zero for all observations and zero for non-observations */
static double stage1_individcov(const int nreta,
								const double reta[static nreta],
								const STAGE1_PARAMS* const params,
								double* icov)
{
	double* testeta = params->testeta;
	let nonzero = params->nonzero;
	let ievaluate_args = &params->ievaluate_args;
	let nrecord = ievaluate_args->nrecord;
	let popparam = &ievaluate_args->popparam;
	let gradient_step = params->stage1->gradient_step;
	double* eval_msec = params->eval_msec;

	assert(nreta > 0);
	assert(nrecord > 0);
	var f_plus_h = mallocvar(double, nrecord);
	var f_minus_h = mallocvar(double, nrecord);
	var yhatvar_plus_h = mallocvar(double, nrecord);
	var yhatvar_minus_h = mallocvar(double, nrecord);
	var J = gsl_matrix_alloc(params->nobs, nreta);
	let nomega = popparam->nomega;
	assert(gradient_step != 0.);

	let advanfuncs = ievaluate_args->advanfuncs;
	let record = ievaluate_args->record;
	let recordinfo = &advanfuncs->recordinfo;

	double xx[OPENPMX_OMEGA_MAX] = { };
	forcount(j, nreta) {
		memcpy(xx, reta, nreta * sizeof(double));
		let v = reta[j];

		/* use the omega diagonal for step size */
		var step = gradient_step;
		let i = nonzero->rowcol[j];
		let omega_var = icov[i * nomega + i];
		if (omega_var != 0.)
			step = sqrt(omega_var);
		step = fmax(step, gradient_step);
		/* TODO: this still needs to be verified that it is optimal compared
		 * to just using gradient_step, or even along the eigenvectors
		 * of icov. In that case we can get the objective function at the
		 * critical points basically for free. */

		/* step eta forward */
		let above = v + step;
		xx[j] = above;
		unreduce_eta(testeta, xx, nonzero);
		struct timespec t3;
		clock_gettime(CLOCK_MONOTONIC, &t3);
		individual_evaluate(ievaluate_args,
							0,				/* dont save imodel */
							0,				/* dont save predictvars */
							0,				/* dont save state */
							f_plus_h, yhatvar_plus_h,	/* we need output */
							0, 0);			/* dont need to calculate objfn */
		timespec_duration(&t3, eval_msec);
		*(params->ineval) += 1;

		/* step eta backward */
		let below = v - step;
		xx[j] = below;
		unreduce_eta(testeta, xx, nonzero);
		clock_gettime(CLOCK_MONOTONIC, &t3);
		individual_evaluate(ievaluate_args,
							0,				/* dont save imodel */
							0,				/* dont save predictvars */
							0,				/* dont save state */
							f_minus_h, yhatvar_minus_h,	/* we need output */
							0, 0);			/* dont need to calculate objfn */
		timespec_duration(&t3, eval_msec);
		*(params->ineval) += 1;

		/* calculate derivatives, scaling by yhatvar */
		var iobs = 0;
		const RECORD* ptr = record;
		forcount(k, nrecord) {
			if (RECORDINFO_EVID(recordinfo, ptr) == 0) {
				let dv = RECORDINFO_DV(recordinfo, ptr);
				let upper = (f_plus_h[k] - dv) / sqrt(yhatvar_plus_h[k]);
				let lower = (f_minus_h[k] - dv) / sqrt(yhatvar_minus_h[k]);
				let deriv = (upper - lower) / (above - below);
				gsl_matrix_set(J, iobs, j, deriv);
				++iobs;
			}
			ptr = RECORDINFO_INDEX(recordinfo, ptr, 1);
		}
	}
	free(f_plus_h);
	free(yhatvar_plus_h);
	free(f_minus_h);
	free(yhatvar_minus_h);

	/* for now accumulate the inverse in reducedcov */
	var reducedcov = gsl_matrix_alloc(nreta, nreta);
	let omegainverse = gsl_matrix_const_view_array(nonzero->inversedata, nreta, nreta);
	gsl_matrix_memcpy(reducedcov, &omegainverse.matrix);

	/* we made J such that tJ*J is tGi*invVi*Gi in Term 5 from Bae and Yim */
	/* multiply and accumulate omega inverse, all in one command */
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., J, J, 1., reducedcov);
	gsl_matrix_free(J);

	/* we have added the contribution of the population omega inverse with
	   the contribution from each individual to get the individual covariance,
	   well, the inverse of it. Now we have to invert to get the actual
	   individual covariance matrix */
	let oldhandler = gsl_set_error_handler_off();
	if (gsl_linalg_cholesky_decomp1(reducedcov) != GSL_SUCCESS) {
		forcount(i, nreta) {
			let val = gsl_matrix_get(reducedcov, i, i);
			gsl_matrix_set(reducedcov, i, i, val + 1e-6); 
		}
		gsl_linalg_cholesky_decomp1(reducedcov);
	}
	gsl_set_error_handler(oldhandler);
	let lndet = matrix_lndet_from_cholesky(reducedcov);
	gsl_linalg_cholesky_invert(reducedcov);

	/* update icov now */
	write_icov_from_reduced(icov, nonzero, reducedcov, nomega);
	gsl_matrix_free(reducedcov);

	return lndet;
}

/* NOTE: this function must be thread safe on the level of an individual */
void stage1_thread(INDIVID* const individ,
				   const ADVANFUNCS* const advanfuncs,
				   const POPMODEL* const popmodel,
				   const NONZERO* const nonzero,
				   const OPTIONS* const options,
				   const SCATTEROPTIONS* const scatteroptions)
{
	/* this function optimizes the full individual eta vector and the
	 * icov_lndet term which are both necessary for the objective function.
	 * In addition, the full (not reduced) individual covariance matrix
	 * is updated */
	let nomega = popmodel->nomega;
	let nreta = nonzero->n;
	double* icov = individ->icov;

	struct timespec t1;
	clock_gettime(CLOCK_MONOTONIC, &t1);

	/* reset result to nothing */
	individ->obs_min2ll = 0.;
	individ->obs_lndet = 0.;
	individ->eta_min2ll = 0.;
	individ->icov_lndet = 0.;
	individ->iobjfn = DBL_MAX;

	/* do the evaluation in a test eta and copy back at the end */
	/* this has to be done because we may unreduce the eta and there we
		assume that it is OPENPMX_OMEGA_MAX wide, whereas in the INDIVID
		it is only nomega wide. We just have to copy back and forth. */
	int stage1_ineval = 0;
	double testeta[OPENPMX_OMEGA_MAX] = { };
	forcount(i, OPENPMX_OMEGA_MAX)
		testeta[i] = NAN;
	let stage1_params = (const STAGE1_PARAMS) {
		.testeta = testeta,
		.nonzero = nonzero,
		.ineval = &stage1_ineval,
		.stage1 = &options->estimate.stage1,
		.ievaluate_args = ievaluate_args_init(individ->record,
											  individ->nrecord,
											  advanfuncs,
											  popmodel->theta,
											  popmodel->ntheta,
											  testeta,
											  popmodel->nomega,
											  popmodel->sigma,
											  popmodel->nsigma,
											  scatteroptions ? scatteroptions->logstream : 0),
		.nobs = individ->nobs,
		.eval_msec = &individ->eval_msec,
	};
	double reta[OPENPMX_OMEGA_MAX] = { };

	/* if no etas or observations, we cant do anything */
	/* we dont have anything to popmodel we just set all etas are zero
	 * because we know that is the minimum 
	 * so it might be interesting to make that possible as well */
	var ieta = individ->eta;
	bool all_eta_zero = true;
	if (nreta == 0 || individ->nobs == 0) {
		memset(ieta, 0, nomega * sizeof(double));
		memset(icov, 0, nomega * nomega * sizeof(double));
		memset(testeta, 0, nomega * sizeof(double)); /* we evaluate at zero eta */

	/* we have etas so we can minimize, start with the current value */
	} else {
		/* start with current value TODO: make this same across runs via besteta */ 
		memcpy(testeta, ieta, nomega * sizeof(double));
		
/// The inner (Stage 1) optimization only optimizes the first, second, 
/// and third terms in the objective function. The fourth term is not
/// dependant on the individual and the fifth term is calculated at the
/// minimum of the first three terms.

		/* optimize the individual eta, the result is written into stage1_params.testeta
		 * (which points to testeta right now) and is also written in reta as well and
		 * write the result back into the individual */
		all_eta_zero = estimate_individual_posthoc_eta(reta, &stage1_params);
		memcpy(ieta, testeta, nomega * sizeof(double));

		/* eta likelihood for part of the objective function */
		/* functions for Bae and Yim objective function Term 3 */
		individ->eta_min2ll = sample_min2ll(nreta, reta, nonzero);
	}

	/* get all the results at the final eta results, fill in the individual yhat and yhatvar */
	/* TODO: maybe dont fill in the values if we will sample the ICOV later?
	   but we do need yhatvar for the calculation of the covariance second term */
	struct timespec t3;
	clock_gettime(CLOCK_MONOTONIC, &t3);
	double obs_min2ll = 0.;
	double obs_lndet = 0.;
	individual_evaluate(&stage1_params.ievaluate_args,
						individ->imodel,
						individ->predictvars,
						individ->istate,
						individ->yhat,
						individ->yhatvar,
						&obs_lndet, &obs_min2ll);
	timespec_duration(&t3, &individ->eval_msec);
	individ->obs_min2ll = obs_min2ll;
	individ->obs_lndet = obs_lndet;
	++stage1_ineval;

	/* we cant do covariance matrix if we have no etas or observations */
	if (nreta == 0 || individ->nobs == 0)
		return;

	/* compute inverse covariance of best fit */
	/* fill in the total individual covariance matrix */
	/* This needs the individual YHATVAR to be calculated whcih we
	 * did in the evaluate above */
	
	/* on first run, do it more than once to stabilize step sizes */
	let niter = all_eta_zero ? 4 : 1;
	forcount(i, niter) {

		/* regular icov calcualtion via first derivatives. This adds the
		 * contribution of the population omega inverse with the
		 * contribution from each individual to get the inverse individual
		 * covariance and then inverts that and writes icov back */
		let last_icov_lndet = individ->icov_lndet;
		individ->icov_lndet = stage1_individcov(nreta,
												reta,
												&stage1_params,
												icov);

		/* if icov does not change, we can stop now */
		if (fabs(last_icov_lndet - individ->icov_lndet) < 0.01)
			break;
	}

	individ->ineval += stage1_ineval;
	timespec_duration(&t1, &individ->stage1_msec);
}

