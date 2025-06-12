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
#include "scatter.h"
#include "advan/advan.h"
#include "utils/c22.h"
#include "utils/various.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

/*--------------------------------------------------------------------*/
/* different optimizers */
#define OPTIMIZER_INNER_BOBYQA

/* different ways to calculate sample min2ll */
//#define SAMPLE_MIN2LL_FROM_INVERSE
#define SAMPLE_MIN2LL_FROM_CHOLESKY

/*--------------------------------------------------------------------*/
static void reduce_eta(double* reta,
					   const double fulleta[static OPENPMX_OMEGA_MAX],
					   const NONZERO* const nonzero)
{
	forcount(i, nonzero->n) {
		let col = nonzero->rowcol[i];
		reta[i] = fulleta[col];
	}
}

static void unreduce_eta(double fulleta[static OPENPMX_OMEGA_MAX],
						 const double* reta,
						 const NONZERO* const nonzero)
{
	memset(fulleta, 0, OPENPMX_OMEGA_MAX * sizeof(double));
	forcount(i, nonzero->n) {
		let col = nonzero->rowcol[i];
		var v = reta[i];
		fulleta[col] = v;
	}
}

#ifdef SAMPLE_MIN2LL_FROM_CHOLESKY
static double sample_min2ll_from_cholesky(const double* const data,
										  const gsl_matrix* const chol)
{
	let n = chol->size1;
	let x = gsl_vector_const_view_array(data, n);

	/* Basically this 'decorrelates' the sample with respect to the covariance matrix
	 (via its cholesky decomposition) which gives the independent samples. The sum of
	 the square is the likelihood */
	double rdata[OPENPMX_OMEGA_MAX];
	var r = gsl_vector_view_array(rdata, n);
	gsl_vector_memcpy(&r.vector, &x.vector);
	gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, chol, &r.vector);

	double sum = 0.;
	forcount(i, n) {
		let v = gsl_vector_get(&r.vector, i);
		sum += v * v;
	}
	return sum;
}
#endif

#ifdef SAMPLE_MIN2LL_FROM_INVERSE
static double sample_min2ll_from_inverse(const double* const data,
										 const gsl_matrix* const inverse)
{
	let n = inverse->size1;
	double ydata[OPENPMX_OMEGA_MAX];
	assert(n < OPENPMX_OMEGA_MAX);
	var y = gsl_vector_view_array(ydata, n);
	let x = gsl_vector_const_view_array(data, n);

	gsl_blas_dgemv(CblasNoTrans, 1., inverse, &x.vector, 0., &y.vector);

//	var lik = 0.;
//	gsl_blas_ddot(&x.vector, &y.vector, &lik);
//	return lik;
	double sum = { };
	forcount(i, n) {
		let x = gsl_vector_get(&x.vector, i);
		let y = gsl_vector_get(&y.vector, i);
		sum += x * y;
	}
	return sum;
}
#endif

/* functions for Bae and Yim objective function Term 3 */
static double sample_min2ll(const int nreta,
							const double reta[static nreta],
							const NONZERO* const nonzero)
{
	assert(nreta == nonzero->n);
	var term3 = 0.;
#ifdef SAMPLE_MIN2LL_FROM_INVERSE
	let omegainverse_data = nonzero->inversedata;
	let omegainverse = gsl_matrix_const_view_array(omegainverse_data, nreta, nreta);
	term3 += sample_min2ll_from_inverse(reta, &omegainverse.matrix);
#endif
#ifdef SAMPLE_MIN2LL_FROM_CHOLESKY
	let omegacholesky_data = nonzero->choleskydata;
	let omegacholesky = gsl_matrix_const_view_array(omegacholesky_data, nreta, nreta);
	term3 += sample_min2ll_from_cholesky(reta, &omegacholesky.matrix);
#endif

#ifdef SAMPLE_MIN2LL_FROM_INVERSE
#ifdef SAMPLE_MIN2LL_FROM_CHOLESKY
	term3 /= 2.;
#endif
#endif
//	assert(isfinite(term3) == 1);
	return term3;
}

typedef struct {
	double* const testeta;
	const NONZERO* const nonzero;
	int* const neval;
	const STAGE1CONFIG* const stage1;
	const IEVALUATE_ARGS ievaluate_args;
	double* const eval_msec;
} STAGE1_PARAMS;

static double stage1_evaluate_individual_iobjfn(const long int nreta,
												const double reta[static nreta],
												void* const data)
{
	let stage1_params = (const STAGE1_PARAMS * const) data;
	let ievaluate_args = &stage1_params->ievaluate_args;

	/* reta and nreta are the variable in the reduced eta/omega matricies */
	/* have to un-reduce the eta values because user code needs the full eta */
	/* the user code can see the testeta because ievaluate_args->popparam->eta
	 * points to the same place */
	unreduce_eta(stage1_params->testeta, reta, stage1_params->nonzero);
	assert(nreta == stage1_params->nonzero->n);

	/* functions for Bae and Yim objective function Term 1 and Term 2 */
	struct timespec t3;
	clock_gettime(CLOCK_REALTIME, &t3);
	let objfn_term1_term2 = individual_fasteval(ievaluate_args);
	timespec_duration(&t3, stage1_params->eval_msec);
	*(stage1_params->neval) += 1;

	/* functions for Bae and Yim objective function Term 3 */
	var term3 = sample_min2ll(nreta, reta, stage1_params->nonzero);

	let iobjfn = objfn_term1_term2 + term3;
	return iobjfn;
}

#ifdef OPTIMIZER_INNER_BOBYQA
#include "bobyqa/bobyqa.h"
#endif

/* NOTE: this function must be thread safe and only touch individual data */
static void estimate_individual_posthoc_eta(double reta[static OPENPMX_OMEGA_MAX],
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
		lower[i] = -DBL_MAX;
		upper[i] = DBL_MAX;
		if (reta[i] != 0.)
			all_eta_zero = false;
	}

	/* extra initial optimization when starting from 0 */
	let n = nreta;
	assert(n < OPENPMX_OMEGA_MAX);
	let neval = stage1->maxeval;

#ifdef OPTIMIZER_INNER_BOBYQA
	int retcode;
	let iprint = 0;
	let npt_recommended = 2*n+1;	/* recommended, used for refine stage */
//	let npt_min = n+2;				/* minimum */
//	let npt_max = (n+1)*(n+2)/2;	/* maximum */

	var npt = npt_recommended;
	let wsize = (npt+5)*(npt+n)+3*n*(n+5)/2 + 10; 	/* a little bit extra room to be sure */
	var w = mallocvar(double, wsize);
	var rhobeg = stage1->step_initial;
	var rhoend = stage1->step_refine;
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
			warning(0, "ID %f eta (initial) not successful (%i)\n", id, retcode);
		}
	}

	/* regular refinement optimization */
	/* we dont need to realloc because the previous memory in w is enough */
	npt = npt_recommended;
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
		warning(0, "ID %f eta (refine) not successful (%i)\n", id, retcode);
	}
	free(w);
#endif

	/* final results, i.e. the etas of the reduced matrix, get expanded */
	assert(stage1_params->testeta == popparam->eta);

	unreduce_eta(stage1_params->testeta, reta, nonzero);
}

/* calculate individual variance covariance matrix
 * yhatvar must be non-zero for all observations and zero for non-observations */
static void stage1_reducedicov(gsl_matrix * const reducedicov,
							   const int nreta,
							   const double reta[static nreta],
							   const STAGE1_PARAMS* const params,
							   const double* icov)
{
	double* testeta = params->testeta;
	let nonzero = params->nonzero;
	let ievaluate_args = &params->ievaluate_args;
	let nrecord = ievaluate_args->nrecord;
	let popparam = &ievaluate_args->popparam;
	let gradient_step = params->stage1->gradient_step;
	double* eval_msec = params->eval_msec;

	let omegainverse_data = nonzero->inversedata;
	let omegainverse = gsl_matrix_const_view_array(omegainverse_data, nreta, nreta);
	gsl_matrix_memcpy(reducedicov, &omegainverse.matrix);

	assert(nreta > 0);
	assert(nrecord > 0);
	var f_plus_h = mallocvar(double, nrecord);
	var f_minus_h = mallocvar(double, nrecord);
	var yhatvar_plus_h = mallocvar(double, nrecord);
	var yhatvar_minus_h = mallocvar(double, nrecord);
	var J = gsl_matrix_alloc(nrecord, nreta);
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
		double step = gradient_step;
		let i = nonzero->rowcol[j];
		let omega_var = icov[i * nomega + i];
		if (omega_var != 0.)
			step = sqrt(omega_var);
		if (step < gradient_step)
			step = gradient_step;
		/* TODO: this still needs to be verified that it is optimal compared to just using gradient_step */

		/* step eta forward */
		let above = v + step;
		xx[j] = above;
		unreduce_eta(testeta, xx, nonzero);
		struct timespec t3;
		clock_gettime(CLOCK_REALTIME, &t3);
		individual_evaluate(ievaluate_args,
							0,				/* dont save imodel */
							0,				/* dont save predictvars */
							0,				/* dont save state */
							f_plus_h, yhatvar_plus_h,	/* we need output */
							0, 0);			/* dont need to calculate objfn */
		timespec_duration(&t3, eval_msec);
		*(params->neval) += 1;

		/* step eta backward */
		let below = v - step;
		xx[j] = below;
		unreduce_eta(testeta, xx, nonzero);
		clock_gettime(CLOCK_REALTIME, &t3);
		individual_evaluate(ievaluate_args,
							0,				/* dont save imodel */
							0,				/* dont save predictvars */
							0,				/* dont save state */
							f_minus_h, yhatvar_minus_h,	/* we need output */
							0, 0);			/* dont need to calculate objfn */
		timespec_duration(&t3, eval_msec);
		*(params->neval) += 1;

		/* calculate derivatives, scaling by yhatvar */
		const RECORD* ptr = record;
		forcount(k, nrecord) {
			let evid = RECORDINFO_EVID(recordinfo, ptr);
			var deriv = 0.;
			if (evid == 0) {
				let dv = RECORDINFO_DV(recordinfo, ptr);
				let upper = (f_plus_h[k] - dv) / sqrt(yhatvar_plus_h[k]);
				let lower = (f_minus_h[k] - dv) / sqrt(yhatvar_minus_h[k]);
				deriv = (upper - lower) / (above - below);
			}
			gsl_matrix_set(J, k, j, deriv);

			ptr = RECORDINFO_INDEX(recordinfo, ptr, 1);
		}
	}
	/* we made J such that tJ*J is tGi*invVi*Gi in Term 5 from Bae and Yim */
	/* multiply and accumulate omega inverse, all in one command */
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., J, J, 1., reducedicov);

	free(f_plus_h);
	free(yhatvar_plus_h);
	free(f_minus_h);
	free(yhatvar_minus_h);
	gsl_matrix_free(J);
}

/* Here we refine the covariance matrix by taking samples at the
 * "sigma" points (1 SD from center) and estimating ln(det(icov))
 * via the Akaike weights. Maybe we could also use the cholesky
 * points, but then the specific points used will depend on the ordering
 * of the etas and that is not desirable for an optimization method */
static double stage1_icov_resample(const gsl_matrix * const reducedicov,
								   double* icovweight,
								   double* icovsample,
								   const int nreta,
								   const double reta[static nreta],
								   const double base_iobjfn,
								   const STAGE1_PARAMS* const stage1_params)
{
	/* eigenvalue decomposition */
	var ework = gsl_matrix_alloc(nreta, nreta);
	var eval = gsl_vector_alloc(nreta);
	var evec = gsl_matrix_alloc(nreta, nreta);
	var w = gsl_eigen_symmv_alloc(nreta);
	var retavals = mallocvar(double, (nreta * 2) * nreta);
	gsl_matrix_memcpy(ework, reducedicov);
	gsl_eigen_symmv(ework, eval, evec, w);

	/* There will be 2*nreta samples and weights and the allocated memory
	 * will allow 2*nomega. The unused ones will be denoted by zero weights
	 * and zero sample etas */
	let nomega = stage1_params->ievaluate_args.popparam.nomega;
	memset(icovweight, 0, 2 * nomega * sizeof(double));
	memset(icovsample, 0, 2 * nomega * nomega * sizeof(double));

	/* make sure we can evaluate at a test point without changing the
	 * predictions of the individual */
	double testreta[OPENPMX_OMEGA_MAX] = { };
	let stage1 = stage1_params->stage1;
	forcount(i, nreta) {
		/* magnitude of step is the square root of eigenvalue */
		var stepsize = sqrt(gsl_vector_get(eval, i));
		if (stepsize < stage1->gradient_step)
			stepsize = stage1->gradient_step;

		/* sample on the positive direction of the eigenvectors and correct for sampling weight */
		forcount(j, nreta) {
			let v = gsl_matrix_get(evec, j, i); /* eigevnectors are coloumns */
			testreta[j] = reta[j] + stepsize * v;
			retavals[(i * 2) * nreta + j] = testreta[j];
		}
		var etaval_iobjfn = stage1_evaluate_individual_iobjfn(nreta, testreta, (void*)stage1_params);
		var delta = etaval_iobjfn - base_iobjfn;
		var lik = exp(-0.5*delta);
		let w1 = lik / exp(-0.5*1.);
			
		/* save the eta used in the individual */
		icovweight[i * 2] = w1;
		forcount(k, nomega)
			icovsample[(i * 2) * nomega + k] = stage1_params->testeta[k];

		/* sample on the negative direction of the eigenvectors and correct for sampling weight */
		forcount(j, nreta) {
			let v = gsl_matrix_get(evec, j, i); /* eigevnectors are coloumns */
			testreta[j] = reta[j] - stepsize * v;
			retavals[(i * 2 + 1) * nreta + j] = testreta[j];
		}
		etaval_iobjfn = stage1_evaluate_individual_iobjfn(nreta, testreta, (void*)stage1_params);
		delta = etaval_iobjfn - base_iobjfn;
		lik = exp(-0.5*delta);
		let w2 = lik / exp(-0.5*1.);

		/* save the eta used in the individual */
		icovweight[i * 2 + 1] = w2;
		forcount(k, nomega)
			icovsample[(i * 2 + 1) * nomega + k] = stage1_params->testeta[k];
	}
	
	/* should we do weighted covariance calculation of the sampled points?
	 * NO because the lndet_icov does not match the sampling matrix!
	 * It should be in the documentation that this is the symmetric matrix
	 * used for sampling and not the one that determined log(det) */
	/* For now we just keep the sampling icov in the same place */
/*	forcount(j, nreta) {
		forcount(k, j + 1) {
			var s = 0.;
			forcount(i, nreta * 2) {
				let xbarj = reta[j];
				let xbark = reta[k];
				let xij = gsl_matrix_get(retavals, i, j);
				let xik = gsl_matrix_get(retavals, i, k);
				let wi = icovweight[i] * 0.5;
				s += wi * (xij - xbarj) * (xik - xbark);
			}
			gsl_matrix_set(reducedicov, j, k, s);
			gsl_matrix_set(reducedicov, k, j, s);
		}
	} */
	
	/* we have to average the weights on the eigenvectors
	 * and take advantage of the fact that determinant of a matrix is
	 * the product of its eigenvectors, in this case, scaled eigenvectors,
	 * with two of them (forward and backward) for each dimension.
	 * Here we take advantage of log(a*b) = log(a) + log(b) which is
	 * probably more accurate. Does this mean that sampling from icov is
	 * not the same as sampling from the individual because our determinant
	 * is asymmetric? I think this does mean the assumption of normal
	 * distribution is relaxed with respect to the minimum objective
	 * function that we will minimize to */
	var sumlogeval = 0.;
	forcount(i, nreta) {
		let w1 = icovweight[i * 2];
		let w2 = icovweight[i * 2 + 1];
		let eivar = gsl_vector_get(eval, i);

		/* weighting via SD */
		if (1) {
			let eisd1 = sqrt(eivar);						/* distance in positive direction */
			let eisd2 = sqrt(eivar); 						/* distance in negative direction */
			let eiwgtsd = (eisd1 * w1 + eisd2 * w2) / 2.;	/* average the weights, i.e, two of them weighted by 0.5 */
			sumlogeval += 2.*log(eiwgtsd);
			
		/* weighting via variance */
		} else {
			let eiwgtvar = (eivar * w1 + eivar * w2) / 2.;	/* average the variances, i.e, two of them weighted by 0.5 */
			sumlogeval += log(eiwgtvar);
		}
	}
	let icov_lndet = -1. * sumlogeval; /* we get lndet from inverse, so we have -1 here */

	free(retavals);
	gsl_eigen_symmv_free(w);
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	gsl_matrix_free(ework);

	return icov_lndet;
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
	clock_gettime(CLOCK_REALTIME, &t1);

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
	let stage1_params = (const STAGE1_PARAMS) {
		.testeta = testeta,
		.nonzero = nonzero,
		.neval = &stage1_ineval,
		.stage1 = &options->estimate.stage1,
		.ievaluate_args = {
			.record = individ->record,
			.nrecord = individ->nrecord,
			.advanfuncs = advanfuncs,
			.popparam = popparam_init(popmodel, advanfuncs, testeta),
			.logstream = scatteroptions->logstream,
		},
		.eval_msec = &individ->eval_msec,
	};
	double reta[OPENPMX_OMEGA_MAX] = { };

	/* if no etas or observations, we cant do anything */
	/* we dont have anything to popmodel we just set all etas are zero
	 * because we know that is the minimum 
	 * so it might be interesting to make that possible as well */
	var ieta = individ->eta;
	if (nreta == 0 || individ->nobs == 0) {
		memset(ieta, 0, nomega * sizeof(double));
		memset(icov, 0, nomega * nomega * sizeof(double));

	/* we have etas so we can minimize, start with the current value */
	} else {
		/* start with current value TODO: make this same across runs via besteta */ 
		memcpy(testeta, ieta, nomega * sizeof(double));

		/* optimize the individual eta, the result is written into stage1_params.testeta
		 * (which points to testeta right now) and is also written in reta as well and
		 * write the result back into the individual */
		estimate_individual_posthoc_eta(reta, &stage1_params);
		memcpy(ieta, testeta, nomega * sizeof(double));

		/* eta likelihood for part of the objective function */
		/* functions for Bae and Yim objective function Term 3 */
		individ->eta_min2ll = sample_min2ll(nreta, reta, nonzero);
	}

	/* get all the results at the final eta results, fill in the individual yhat and yhatvar */
	/* TODO: maybe dont fill in the values if we will sample the ICOV later?
	   but we do need yhatvar for the calculation of the covariance second term */
	struct timespec t3;
	clock_gettime(CLOCK_REALTIME, &t3);
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
	/* fill in the total individual covariance matrix from the reduced one */
	/* This needs the individual YHATVAR to be calculated */
	var reducedicov = gsl_matrix_alloc(nreta, nreta);
	stage1_reducedicov(reducedicov,
					   nreta,
					   reta,
					   &stage1_params,
					   icov);
	/* we have added the contribution of the population omega inverse with
	   the contribution from each individual to get the individual covariance,
	   well, the inverse of it. Now we have to invert to get the actual
	   individual covariance matrix */
	/* This could possibly be safer with another decomposition method?
	   maybe it does not matter */
	gsl_linalg_cholesky_decomp1(reducedicov);
	individ->icov_lndet = matrix_lndet_from_cholesky(reducedicov);

	if (!isfinite(individ->icov_lndet)) {
		warning(0, "individual lndet is not finite\n"); /* TODO: add info about ID and eta */
		individ->icov_lndet = 100.;
	}
//	assert(isfinite(individ->icov_lndet) == 1);
	gsl_linalg_cholesky_invert(reducedicov);

	/* Here we refine the covariance matrix by taking samples at the
	 * inflection points (1 SD from center) and doing a weighted covariance
	 * matrix estimation. */
	if (options->estimate.stage1.icov_resample) {
		let base_iobjfn = obs_min2ll + obs_lndet + individ->eta_min2ll;
		let icov_resample_lndet = stage1_icov_resample(reducedicov,
													   individ->icovweight,
													   individ->icovsample,
													   nreta,
													   reta,
													   base_iobjfn,
													   &stage1_params);
		if (isfinite(icov_resample_lndet))
			individ->icov_lndet = icov_resample_lndet;
		else
			warning(0, "icov resample lndet is not finite, ignoring\n");
	}

	/* We use icov for determining the step size when doing the icov calculation
	   via gradients in stage1_reducedicov so we cant skip this */
	let rowcol = nonzero->rowcol;
	forcount(i, nreta) {
		forcount(j, nreta) {
			let row = rowcol[i];
			let col = rowcol[j];
			let v = gsl_matrix_get(reducedicov, i, j);
			icov[row * nomega + col] = v;
		}
	}

	gsl_matrix_free(reducedicov);

	individ->ineval += stage1_ineval;
	timespec_duration(&t1, &individ->stage1_msec);
}

