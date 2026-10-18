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
 
/// This file calculates the covariance matrix of the estimate to provide
/// information about upper and lower bounds of paramater estimates.
///
/// For off-diagonal values it is not really clear how to transform these from the
/// unbounded encoded space to the bounded paramater space, so no limits
/// are calculated for them.
///
/// All work is done in the encoded (transformed) parameter space using
/// central finite differences of the objective function (-2LL).
/// Two matrices are relevant:
///
/// - **R**: the Hessian of the total objfn (NONMEM's R matrix).
/// - **S**: the OPG/BHHH cross-product of per-individual gradients (S = J'J).
///
/// They are combined into the sandwich estimator `Cov = Rinv * S * Rinv`, or
/// used alone (`Sinv` or `Rinv`, scaled because the objfn is -2LL). 
/// Quadratic smoothing can then refine the result. The covariance type is
/// chosen with `COVARIANCE_S_MATRIX`, `COVARIANCE_R_MATRIX`,
/// `COVARIANCE_SANDWICH` or the default `COVARIANCE_SANDWICH_QUADRATIC`.
///
///   Results go to `<filename>.cov` and `<filename>.cov.ext` when the
///   source has a filename.
/// - `pmx_covariance_limit(pmx, param)`: looks up the stored `COVLIMIT`
///   for a given parameter type and index. If no covariance has
///   been calculated an empty `COVLIMIT` is returned.
///
/// ## Reporting
///
/// - Correlation pairs sorted by absolute value are reported.
/// - Eigen-analysis of S and/or R finds all parameters confined to the
///   near-null subspace (eigenvalues below `reltol * lambda_max`), also
///   reporting negative eigenvalues.
/// - When diagnostic samples are requested the predicted vs observed
///   dobjfn is written to the diagnostics file and summarised
///   (log-ratio mean and sd).

#include "openpmx.h"
#include "encode.h"
#include "print.h"
#include "scatter.h"
#include "stage1.h"
#include "pmxstate.h"
#include "defines.h"
#include "githash.h"
#include "linalg.h"
#include "advan/advan.h"
#include "utils/c22.h"
#include "utils/vector.h"
#include "utils/various.h"

#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_eigen.h>

#define S_FLOOR_RELTOL 1e-7
#define R_FLOOR_RELTOL 1e-5

#define COVARIANCE_DOBJFN_TARGET_LOWER		4.
#define COVARIANCE_DOBJFN_TARGET_UPPER		8.

typedef struct {
	const double* offset;
	double objfn;
} VB_SAMPLE;

typedef struct {
	VECTOR(VB_SAMPLE) samples;
	const int nparam;
} VB_HISTORY;

static void history_append(VB_HISTORY* const hist,
						   const double* offset,
						   const double objfn)
{
	var p = mallocvar(double, hist->nparam);
	memcpy(p, offset, hist->nparam * sizeof(double));
	
	vector_append(hist->samples, 
		(VB_SAMPLE) {
			.offset = p,
			.objfn = objfn,
		});
}

static VB_HISTORY history_new(const int nparam)
{
	return (VB_HISTORY) {
		.samples = { { 0 } },
		.nparam = nparam,
	};
}	

static void history_free(VB_HISTORY* const hist)
{
	forvector_ptr(p, hist->samples)
		free((void*)p->offset);
	vector_free(hist->samples);
}

static int history_prune_above(VB_HISTORY* const hist, const double threshold, FILE* outstream)
{
	int kept = 0;
	let oldsize = hist->samples.size;
	forcount(s, oldsize) {
		let sp = &hist->samples.ptr[s];
		if (!(sp->objfn <= threshold)) { /* handles NAN peroperly */
			info(outstream, "prune %f\n", sp->objfn);
			free((void*)sp->offset);
		} else {
			if (kept != s)
				hist->samples.mutptr[kept] = *sp;
			kept++;
		}
	}
	vector_resize(hist->samples, kept);
	return oldsize - kept;
}

/*
 * Fit a quadratic surface with fixed minimum f0 at center and offset 0:
 *   f(x) - f0 = (1/2) x^T H x
 */
static void vb_fit_quadratic(const VB_HISTORY* const hist,
							 const double f0,
							 FILE* outstream,
							 gsl_matrix* const cov)
{
	let ndim = hist->nparam;
	let nsamples = hist->samples.size;

	/* Number of quadratic basis terms: ndim*(ndim+1)/2 (no intercept, no linear terms) */
	let nfitparam = ndim * (ndim + 1) / 2;
	if (nsamples <= nfitparam) {
		info(outstream, "failed need %i samples, have %i\n", (int)nfitparam, (int)nsamples);
		return;
	}

	/* Build design matrix X [nsamples x nfitparam] and response vector y */
	int rc = 0;
	double chisq = -1.;
	var c = gsl_vector_alloc(nfitparam);
	{
		gsl_matrix* X = gsl_matrix_alloc(nsamples, nfitparam);
		gsl_vector* y = gsl_vector_alloc(nsamples);
		{
			var x = mallocvar(double, ndim);
			forcount(s, nsamples) {
				/* centered coordinates */
				let sp = &hist->samples.ptr[s];
				forcount(j, ndim)
					x[j] = sp->offset[j];

				/* quadratic terms only */
				int col = 0;
				forcount(i, ndim) {
					for (int j = i; j < ndim; j++) {
						let val = (i == j) ? 0.5 * x[i] * x[j] : x[i] * x[j];
						gsl_matrix_set(X, s, col++, val);
					}
				}
				assert(col == nfitparam);

				/* Shift response vector by f0 */
				/* f0 is assumed to be the true objfn at offset 0. If it is the
				 * estimation-time value rather than one evaluated the same way
				 * as the samples, the fit has a systematic offset. */
				gsl_vector_set(y, s, sp->objfn - f0);
			}
			free(x);
		}
		/* least squares without intercept via GSL multifit */
		{
			var wcov  = gsl_matrix_alloc(nfitparam, nfitparam);
			var work = gsl_multifit_linear_alloc(nsamples, nfitparam);
			/* NOTE: unweighted. Samples far from the optimum (up to f0+30
			 * after pruning) have the largest |residuals| if the surface is
			 * not exactly quadratic, so they dominate the fit although the
			 * Hessian at the optimum is what is wanted. gsl_multifit_wlinear with e.g.
			 * w = 1/(1 + objfn - f0) would favour the near points but this
			 * requires more samples for the same amount of information. */
			rc = gsl_multifit_linear(X, y, c, wcov, &chisq, work);
			gsl_multifit_linear_free(work);
			gsl_matrix_free(wcov);
		}
		gsl_matrix_free(X);
		gsl_vector_free(y);
	}
	if (rc != 0) {
		info(outstream, "failed (%i)\n", rc);
		gsl_vector_free(c);
		return;
	}

	info(outstream, "nsamples %i chisq/(nsamples-nfitparam) %f\n",
		nsamples,
		chisq/(nsamples-nfitparam));
	
	/* extract Hessian H directly from fitted coefficients c */
	var hessian = gsl_matrix_alloc(ndim, ndim);
	gsl_matrix_set_zero(hessian);
	var col = 0;
	forcount(i, ndim) {
		for (int j = i; j < ndim; j++) {
			let v = gsl_vector_get(c, col++);
			if (i == j) {
				gsl_matrix_set(hessian, i, i, v);
			} else {
				gsl_matrix_set(hessian, i, j, v);
				gsl_matrix_set(hessian, j, i, v);
			}
		}
	}
	gsl_vector_free(c);

	/* invert hessian and scale in-place before copying to user */
	var saved_hessian = gsl_matrix_alloc(ndim, ndim);
	gsl_matrix_memcpy(saved_hessian, hessian);
	if (gsl_linalg_cholesky_decomp1(hessian) != GSL_SUCCESS) {
		info(outstream, "cov quadratic hessian not positive-definite, failing\n");
		goto done;
	}
	if (gsl_linalg_cholesky_invert(hessian) != GSL_SUCCESS) {
		info(outstream, "cov quadratic cholesky invert failed\n");
		goto done;
	}

	/* success, so copy results to user */
	/* hessian here is the inverse Hessian of the -2LL surface, and
	 * Cov = 2 * H^-1 on that scale (f - f0 = x' Cov^-1 x = 1/2 x'Hx). */
	gsl_matrix_scale(hessian, 2.);
	gsl_matrix_memcpy(cov, hessian);
	
done:
	gsl_matrix_free(saved_hessian);
	gsl_matrix_free(hessian);
}

#define HESSIAN_MSTEP_MAXITER				5

static double encode_covariance_evaluate(ENCODE* const test,
										 IDATA* const idata,
										 const ADVANFUNCS* const advanfuncs,
										 const OPTIONS* const options)
{
	var popmodel = &test->popmodel;
	let omegainfo = &test->omegainfo;
	let nonzero = &omegainfo->nonzero;
	var scatteroptions = (SCATTEROPTIONS) {
		.stage1_order = true,
		.stage1_evaluate = true,
	};
	scatter_threads(idata, advanfuncs, popmodel, nonzero, options, &scatteroptions, stage1_thread);

	let objfn = idata_objfn(idata, omegainfo->omega_nonzero_lndet);
	popmodel->result.objfn = objfn;
	popmodel->result.type = OBJFN_CURRENT;
	popmodel->result.neval += 1;
	popmodel->result.nsig = 0.;

	return objfn;
}

static double encode_covariance_evaluate_offset(ENCODE* const encode,
												const double* offset,
												IDATA* const idata,
												const IDATAETAS* const warm_etas,
												const ADVANFUNCS* const advanfuncs,
												VB_HISTORY* history,
												const OPTIONS* const options)
{
	idata_etas_write(idata, warm_etas);
	encode_untransform(encode, offset);
	let objfn = encode_covariance_evaluate(encode, idata, advanfuncs, options);
	if (history)
		history_append(history, offset, objfn);
		
	return objfn;
}

static void show_progress(const POPMODEL* const popmodel,
						  const struct timespec* const begin,
						  const bool details,
						  FILE* outstream,
						  FILE* extstream,
						  const char* extra)
{
	struct timespec now;
	clock_gettime(CLOCK_MONOTONIC, &now);
	let t = timespec_time_difference(begin, &now) / 1000.;
	popmodel_eval_information(popmodel, t, 0, details, outstream, extstream, extra);
}

typedef struct {
	ENCODE* const encode;
	VB_HISTORY* history;
	IDATA* const idata;
	const ADVANFUNCS* const advanfuncs;
	const OPTIONS* const options;
	const bool details;
	const IDATAETAS* const warm_etas;
	const double f0;
	const struct timespec* const begin;
	FILE* const outstream;
	FILE* const extstream;
} COVARGS;

/* easier to declare an array that can take an offset for ENCODE */
typedef typeof(((ENCODE*)0)->offset) OFFSET;

/* Evaluate the total objfn at scale * direction in the encoded space,
 * starting from the warm etas. The point is recorded in the history
 * (if any) and reported as progress. If iobjfn_out is given it receives
 * the per-individual objfn at this point. */
static double hessian_eval_step(const COVARGS* const covargs,
								const double* const direction,
								const double scale,
								double* const iobjfn_out,
								const char* const message)
{
	var encode = covargs->encode;
	var idata = covargs->idata;

	OFFSET offset = { 0 };
	forcount(k, encode->nparam)
		offset[k] = scale * direction[k];
	let f = encode_covariance_evaluate_offset(encode, 
											  offset, 
											  idata, 
											  covargs->warm_etas, 
											  covargs->advanfuncs, 
											  covargs->history, 
											  covargs->options);
	show_progress(&encode->popmodel,
				  covargs->begin,
				  covargs->details,
				  covargs->outstream,
				  covargs->extstream,
				  message);

	if (iobjfn_out)
		forcount(k, idata->nindivid)
			iobjfn_out[k] = idata->individ[k].iobjfn;
	return f;
}

/* Next step to try, given that step h gave an average objfn change of
 * avgdiff. dobjfn scales as ~0.5*kappa*h^2 locally, so rather than
 * bisecting, solve directly for the h expected to hit the target.
 * Converges in far fewer evaluations than doubling/bisection because it
 * uses the known h^2 scaling instead of just monotonicity.
 *
 * h_low/h_high are kept only as a fallback bracket for cases where the
 * local quadratic approximation fails (e.g. BLQ discontinuities making
 * dobjfn(h) non-monotonic). */
static double hessian_next_step(const double h,
								const double avgdiff,
								const double target_delta_lower,
								const double target_delta_upper,
								double* const h_low,
								double* const h_high)
{
	let target = sqrt(target_delta_lower * target_delta_upper);
	/* NaN means the model could not be evaluated at this step (too
	 * large), so shrink. Exactly zero change means the step was too
	 * small, so grow. */
	var ratio = 4.;
	if (isnan(avgdiff))
		ratio = 0.25;
	else if (avgdiff > 0.)
		ratio = sqrt(target / avgdiff);
	ratio = fmin(fmax(ratio, 0.1), 10.); /* guard a wild jump off one noisy sample */
	var hnext = h * ratio;

	if (avgdiff <= target_delta_lower)
		*h_low = h;
	else
		*h_high = h; /* includes NaN, since the comparison above is false */
	if (*h_low > 0. && *h_high > 0. && (hnext <= *h_low || hnext >= *h_high))
		hnext = (*h_low + *h_high) / 2.; /* jump left the bracket, fall back to bisection */
	return hnext;
}

/* Perturbs along an arbitrary direction. For the S matrix the 
 * canonical-axis Jacobian is used but for the R matrix the steps are 
 * along the eigenvectors of S.
 *
 * Starting from step h, evaluates the objfn at +h and -h along direction
 * and adapts h until the average |objfn - f0| falls inside
 * (target_delta_lower, target_delta_upper). Returns the step used. The
 * outputs (per-individual objfn and total objfn at +h and -h, each
 * optional) always belong to the returned step. */
static double pmx_hessian_stepsize_dir(const COVARGS* const covargs,
                                       const double* const direction,
                                       double h,
                                       double* const iobjfn_plus_out,
                                       double* const iobjfn_minus_out,
                                       double* const fplus_out,
                                       double* const fminus_out,
                                       const double target_delta_lower,
                                       const double target_delta_upper,
                                       const char* message)
{
	let f0 = covargs->f0;
	var h_low = 0.;
	var h_high = 0.;
	var h_used = h;
	forcount(iter, HESSIAN_MSTEP_MAXITER) {
		h_used = h;

		char message2[256];
		snprintf(message2, sizeof(message2), "%s step %g", message, h_used);
		let fplus = hessian_eval_step(covargs, direction, h_used, iobjfn_plus_out, message2);
		let fminus = hessian_eval_step(covargs, direction, -h_used, iobjfn_minus_out, 0);
		if (fplus_out)
			*fplus_out = fplus;
		if (fminus_out)
			*fminus_out = fminus;

		/* done if the average change falls within the target bracket */
		let avgdiff = (fabs(fplus - f0) + fabs(fminus - f0)) / 2.;
		if (avgdiff > target_delta_lower &&
			avgdiff < target_delta_upper)
			return h_used;

		h = hessian_next_step(h_used, avgdiff,
							  target_delta_lower, target_delta_upper,
							  &h_low, &h_high);
	}

	info(covargs->outstream, "warning: eigenmode step size did not reach objfn "
		"change of [%g,%g] after %i iterations, using h=%g\n",
		target_delta_lower, target_delta_upper,
		HESSIAN_MSTEP_MAXITER, h_used);

	return h_used;
}

/* S matrix: the OPG/BHHH cross-product-of-gradients estimate,
 * S = J'J, where J is the per-individual Jacobian of iobjfn with
 * respect to each transformed parameter (central differences).
 * S is kept uninverted here, same as pmx_hessian_R above, so the caller
 * decides whether to invert alone or combine into a sandwich
 * estimator with R.
 *
 * Its difficult to choose the correct dobjfn for the S matrix because
 * it concerns the total objfn but the derivatives are for each individual
 * which might result in too small changes for each individual. The way
 * the scatter works all individuals get the same popmodel so its
 * difficult to imagine an elegant way to make sure the dobjfn is
 * within some range for each individual. */ 
static gsl_matrix* pmx_opg_S(const COVARGS* const covargs,  
							 const ENCODELABEL* labels,
							 const double step_size_init)
{
	var encode = covargs->encode;
	var idata = covargs->idata;
	let nparam = encode->nparam;
	let nindivid = idata->nindivid;

	var J = gsl_matrix_alloc(nindivid, nparam);
	var iobjfn_plus = mallocvar(double, nindivid);
	var iobjfn_minus = mallocvar(double, nindivid);

	/* Compute OPG Jacobian along canonical coordinate axes also known
	 * as the BHHH or OPG/Fisher information estimate */
	forcount(j, nparam) {
		OFFSET direction = { 0 };
		direction[j] = 1.;
		
		char message[128];
		snprintf(message, sizeof(message), " %s", labels[j]);

		let h = pmx_hessian_stepsize_dir(covargs, direction, 
										 step_size_init, 
										 iobjfn_plus, iobjfn_minus,
										 0, 0,
										 COVARIANCE_DOBJFN_TARGET_LOWER,
										 COVARIANCE_DOBJFN_TARGET_UPPER,
										 message);
	
		/* Compute central difference gradient for column j */
		forcount(k, nindivid) {
			let v = (iobjfn_plus[k] - iobjfn_minus[k]) / (2. * h);
			gsl_matrix_set(J, k, j, v);
		}
	}
	free(iobjfn_plus);
	free(iobjfn_minus);

	/* Form S = J^T * J -- BHHH / outer-product-of-gradients approximation */
	var S = gsl_matrix_alloc(nparam, nparam);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., J, J, 0., S);
	gsl_matrix_free(J);

	return S;
}

static double pmx_hessian_R_crossterm(const COVARGS* const covargs,
									  const double* offset,
									  const char* message)
{
	let val = encode_covariance_evaluate_offset(covargs->encode, 
												offset, 
												covargs->idata, 
												covargs->warm_etas, 
												covargs->advanfuncs, 
												covargs->history,
												covargs->options);
	show_progress(&covargs->encode->popmodel, 
				  covargs->begin, 
				  covargs->details, 
				  covargs->outstream, 
				  covargs->extstream, 
				  message);
	return val;
}

static void triple_product(const gsl_matrix* A,
						   const gsl_matrix* B,
						   gsl_matrix* out)
{
	let nparam = (int)out->size1;
	var tmp = gsl_matrix_alloc(nparam, nparam);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., A, B, 0., tmp);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans,   1., tmp, A, 0., out);
	gsl_matrix_free(tmp);
}

/* R matrix: the Hessian of the *total* objective function (summed
 * over all individuals), via central second differences in the
 * transformed parameter space. This is distinct from the OPG/BHHH
 * "S" matrix (S/cov below, built from per-individual score outer
 * products) -- R uses second derivatives of the objfn directly, no
 * per-individual decomposition needed. Combining both gives the
 * robust sandwich estimator Cov = Rinv * S * Rinv. */
static gsl_matrix* pmx_hessian_R(const COVARGS* const covargs,
								 const gsl_matrix* const S,
								 const double step_size_init)
{
	var encode = covargs->encode;
	var idata = covargs->idata;
	var advanfuncs = covargs->advanfuncs;
	var warm_etas = covargs->warm_etas;
	var f0 = covargs->f0;

	let nparam = encode->nparam;
	let _off = advanfuncs->recordinfo._offset1 ? 1 : 0;

	/* rotate: eigendecompose S to get the directions to step along.
	 * Scopy protects the caller's S from gsl_eigen_symmv's in-place
	 * destruction -- S is still needed untouched for the sandwich
	 * estimator after this returns. */
	var evalS = gsl_vector_alloc(nparam);
	var U = gsl_matrix_alloc(nparam, nparam);
	{
		/* if S is present use it, otherwise just guess */
		var Scopy = gsl_matrix_alloc(nparam, nparam);
		if (S) { 
			gsl_matrix_memcpy(Scopy, S);
		} else {
			gsl_matrix_set_identity(Scopy);
			gsl_matrix_scale(Scopy, 1. / (step_size_init * step_size_init)); /* curvature, not step */
		}
		
		var w = gsl_eigen_symmv_alloc(nparam);
		gsl_eigen_symmv(Scopy, evalS, U, w);
		gsl_eigen_symmv_free(w);
		gsl_matrix_free(Scopy);
	}
	let lambda_max = gsl_vector_max(evalS);
	var Rp = gsl_matrix_alloc(nparam, nparam);

	/* diagonal: fresh adaptive search along each eigenvector, seeded
	 * from S eigenvalue via h ~ 2/sqrt(curvature), curvature
	 * ~ eigenvalue/2 per E[R]=E[S]/2 */
	var step_size = mallocvar(double, nparam); /* per eigenmode, filled in when doing diagonal */
	let avg_delta = (COVARIANCE_DOBJFN_TARGET_LOWER + COVARIANCE_DOBJFN_TARGET_UPPER) / 2.;
	{
		double fplus, fminus;
		OFFSET direction;
		forcount(p, nparam) {
			forcount(k, nparam) 
				direction[k] = gsl_matrix_get(U, k, p);

			char message[128];
			snprintf(message, sizeof(message), " R(%i,%i)", p + _off, p + _off);

			let lam = fmax(gsl_vector_get(evalS, p), R_FLOOR_RELTOL * lambda_max);
			let h_init = 2. * sqrt(avg_delta) / sqrt(lam);
			step_size[p] = pmx_hessian_stepsize_dir(covargs, direction, h_init,
													0, 0, 
													&fplus, &fminus,
													COVARIANCE_DOBJFN_TARGET_LOWER, 
													COVARIANCE_DOBJFN_TARGET_UPPER,
													message);
			let h = step_size[p];
			let hpp = (fplus - 2. * f0 + fminus) / (h * h);
			gsl_matrix_set(Rp, p, p, hpp);
		}
	}
	gsl_vector_free(evalS);

	/* cross terms: perturb along hp*v_p + hq*v_q, same 4-evaluation
	 * pattern as canonical-axis cross terms */
	forcount(p, nparam) {
		for (int q = p + 1; q < nparam; ++q) {
			let hp = step_size[p];
			let hq = step_size[q];
			
			char message[128];
			snprintf(message, sizeof(message), " R(%i,%i)", p + _off, q + _off);
		
			/* both positive */
			OFFSET offset = { 0 };
			forcount(k, nparam)
				offset[k] = hp * gsl_matrix_get(U, k, p) + hq * gsl_matrix_get(U, k, q);
			let fpp = pmx_hessian_R_crossterm(covargs, offset, message);

			/* one positive, one negative */
			forcount(k, nparam)
				offset[k] = hp * gsl_matrix_get(U, k, p) - hq * gsl_matrix_get(U, k, q);
			let fpm = pmx_hessian_R_crossterm(covargs, offset, 0);

			/* one negative, one positive */
			forcount(k, nparam)
				offset[k] = -hp * gsl_matrix_get(U, k, p) + hq * gsl_matrix_get(U, k, q);
			let fmp = pmx_hessian_R_crossterm(covargs, offset, 0);

			/* both negative */
			forcount(k, nparam)
				offset[k] = -hp * gsl_matrix_get(U, k, p) - hq * gsl_matrix_get(U, k, q);
			let fmm = pmx_hessian_R_crossterm(covargs, offset, 0);

			/* set off-diagonal of symmetric R matrix based on samples */
			let hpq = (fpp - fpm - fmp + fmm) / (4. * hp * hq);
			gsl_matrix_set(Rp, p, q, hpq);
			gsl_matrix_set(Rp, q, p, hpq);
		}
	}
	free(step_size);

	/* reset the old eta state */
	idata_etas_write(idata, warm_etas);

	/* de-rotate: R = U * Rp * U' back to canonical parameter basis --
	 * from here on R is an ordinary Hessian in the same basis as S,
	 * usable by the existing regularize/sandwich/report code unchanged */
	var R = gsl_matrix_alloc(nparam, nparam);
	triple_product(U, Rp, R);
	gsl_matrix_free(Rp);
	gsl_matrix_free(U);

	return R;
}

static FILE* estimate_covariance_results_fopen(const char* name,
									const char* ext)
{
	char fname[PATH_MAX];
	snprintf(fname, sizeof(fname), "%s%s", name, ext);
	return fopen(fname, "w");
}

static void decode_popmodel_limit(ENCODE* const encode, const gsl_matrix* const cov, const int j, const double z) 
{
	let variance = gsl_matrix_get(cov, j, j);
	let SD = sqrt(variance);

	OFFSET offset = { 0 };
	offset[j] = z * SD;
	
	encode_untransform(encode, offset);
}

static bool cov_from_inverse(const gsl_matrix* const M,
							 const double factor,
							 const char* const what,
							 FILE* outstream,
							 gsl_matrix* const cov)
{
	var tmp = gsl_matrix_alloc(M->size1, M->size2);
	gsl_matrix_memcpy(tmp, M);
	let ok = gsl_linalg_cholesky_decomp1(tmp) == GSL_SUCCESS &&
			 gsl_linalg_cholesky_invert(tmp) == GSL_SUCCESS;
	if (ok) {
		gsl_matrix_scale(tmp, factor);
		gsl_matrix_memcpy(cov, tmp);
	} else {
		info(outstream, "error: singular or ill-conditioned %s matrix\n", what);
	}
	gsl_matrix_free(tmp);
	return ok;
}

static bool covariance_calculate(const gsl_matrix* const R,
								 const gsl_matrix* const S,
								 FILE* outstream,
								 gsl_matrix* const cov)
{
	let nparam = cov->size1;

	/* sandwich estimator: cov = Rinv * S * Rinv */
	if (R && S) {
		info(outstream, "using sandwich\n");
		var Rinv = gsl_matrix_alloc(nparam, nparam);
		gsl_matrix_memcpy(Rinv, R);
		if (gsl_linalg_cholesky_decomp1(Rinv) != GSL_SUCCESS ||
			gsl_linalg_cholesky_invert(Rinv) != GSL_SUCCESS) {
			info(outstream, "warning: singular or ill-conditioned Hessian (R) matrix\n");
			gsl_matrix_free(Rinv);
			return false;
		}
		triple_product(Rinv, S, cov);
		gsl_matrix_free(Rinv);
		return true;
	}
	/* objfn is -2LL: E[R] = 2*I and E[S] = 4*I (so E[R] = E[S]/2),
	 * hence cov = I^-1 = 2*Rinv or 4*Sinv */
	/* Using 2 because it seems to work better in most (but certainly)
	 * not all cases */
	if (S) {
		info(outstream, "using S matrix\n");
		return cov_from_inverse(S, 2., "OPG/BHHH", outstream, cov);
	}
	if (R) {
		info(outstream, "using R matrix\n");
		return cov_from_inverse(R, 2., "R", outstream, cov);
	}
	return false;
}

typedef struct { 
	const char* i;
	const char* j;
	double val;
} CORRENTRY;

static int correntry_sort_callback(const void* _v1, const void* _v2)
{
	let v1 = (const CORRENTRY*)_v1;
	let v2 = (const CORRENTRY*)_v2;
	let fv1 = fabs(v1->val);
	let fv2 = fabs(v2->val);
	if (fv1 < fv2)
		return 1;
	if (fv1 > fv2)
		return -1;
	return 0;
}

static void print_correlations(const gsl_matrix* const corr, const ENCODELABEL* labels, FILE* outstream)
{
	let nparam = corr->size1;

	let entries = mallocvar(CORRENTRY, nparam * (nparam + 1) / 2);
	int k = 0;
	forcount(i, nparam) {
		forcount(j, i) {
			entries[k++] = (CORRENTRY) { .i = labels[i],
										 .j = labels[j],
										 .val = gsl_matrix_get(corr, i, j) };
		}
	}
	qsort(entries, k, sizeof(entries[0]), correntry_sort_callback);

	info(outstream, "correlations\n");
	info(outstream, OPENPMX_SFORMAT OPENPMX_SFORMAT OPENPMX_SFORMAT "\n", "CORRELATION", "PARAMETER1", "PARAMETER2");
	forcount(i, k) {
		let p = &entries[i];
		info(outstream, OPENPMX_FFORMAT OPENPMX_SFORMAT OPENPMX_SFORMAT "\n", p->val, p->i, p->j);
	}
	free(entries);
}

typedef struct {
	const char* name;
	double contrib;
} EIGCONTRIB;

static int eigcontrib_sort_callback(const void* _v1, const void* _v2)
{
	let v1 = (const EIGCONTRIB*)_v1;
	let v2 = (const EIGCONTRIB*)_v2;
	if (v1->contrib < v2->contrib)
		return 1;
	if (v1->contrib > v2->contrib)
		return -1;
	return 0;
}

/* Finds ALL non-identifiable parameters from a raw (unnormalized)
 * symmetric matrix M (S or R) -- not just the single smallest
 * eigenvalue. Sums each parameter's squared component across every
 * eigenvalue <= reltol * lambda_max, i.e. the diagonal of the
 * projector onto that whole near-null subspace. This is
 * basis-independent: within a degenerate or near-degenerate cluster of
 * small eigenvalues, individual eigenvectors are only defined up to
 * rotation, but the projector onto the subspace they jointly span is
 * not -- so this correctly finds multiple simultaneous
 * non-identifiabilities (a fully unused parameter, two parameters only
 * identifiable as a sum, etc.) without depending on which particular
 * basis GSL happened to return. Score is in [0,1]; near 1 means "this
 * parameter's information is (almost) entirely confined to the
 * rank-deficient subspace".
 *
 * Works for both S (OPG, PSD by construction, rank <= N individuals)
 * and R (Hessian, may be indefinite). Negative eigenvalues are counted
 * as null here, but also reported separately, since for R they mean
 * saddle/negative curvature or finite-difference noise rather than a
 * flat direction. */
static void print_nonidentifiable_parameters(const gsl_matrix* const M,
											 const char* const name,
											 const ENCODELABEL* labels,
											 const double reltol,
											 FILE* outstream)
{
	let nparam = (int)M->size1;
	if (nparam == 0)
		return;

	var Mcopy = gsl_matrix_alloc(nparam, nparam);
	gsl_matrix_memcpy(Mcopy, M);
	var eval = gsl_vector_alloc(nparam);
	var evec = gsl_matrix_alloc(nparam, nparam);
	var w    = gsl_eigen_symmv_alloc(nparam);
	gsl_eigen_symmv(Mcopy, eval, evec, w);
	gsl_eigen_symmv_free(w);
	gsl_matrix_free(Mcopy);

	double emin, emax;
	gsl_vector_minmax(eval, &emin, &emax);
	let scale = fmax(fabs(emin), fabs(emax));
	let floor = reltol * scale;

	var entries = mallocvar(EIGCONTRIB, nparam);
	forcount(i, nparam) {
		let v = (EIGCONTRIB) { .name = labels[i],
							   .contrib = 0. };
		entries[i] = v;
	}
	
	int nnull = 0;
	int nneg = 0;
	forcount(k, nparam) {
		let ev = gsl_vector_get(eval, k);
		if (ev < 0.) {
			++nneg;
			continue;
		}
		if (ev > floor)
			continue;
		++nnull;
		forcount(i, nparam) {
			let v = gsl_matrix_get(evec, i, k);
			entries[i].contrib += v * v;
		}
	}
	gsl_matrix_free(evec);
	gsl_vector_free(eval);

	qsort(entries, nparam, sizeof(entries[0]), eigcontrib_sort_callback);

	if (nneg)
		info(outstream, "%s %i of %i eigenvalue(s) negative: not a minimum or noisy Hessian\n",
			 name, nneg, nparam);

	if (nnull) {
		info(outstream, "%s near-null subspace: %i of %i eigenvalue(s) in [0, %g x largest |eigenvalue|]\n",
			 name, nnull, nparam, reltol);

		info(outstream, OPENPMX_SFORMAT OPENPMX_SFORMAT "\n", "PARAMETER", "NULLSPACE_CONTRIB");
		forcount(i, nparam) 
			info(outstream, OPENPMX_SFORMAT OPENPMX_FFORMAT "\n", entries[i].name, entries[i].contrib);
	}

	free(entries);
}

/*
 * S0           : original covariance (n x n, SPD)
 * x            : sample
 * observed_2ll : observed Gaussian -2LL of x
 * S            : output, adjusted covariance with 2LL(S; x) = observed_2ll
 * kl           : output, KL(N(0,S) || N(0,S0)), may be NULL
 *
 * Returns GSL_SUCCESS, GSL_EDOM (S0 not SPD or x == 0)
 */
static int covariance_match_2ll(const gsl_matrix *S0,
								const gsl_vector *x,
								double observed_2ll,   /* target x'S^-1 x */
								gsl_matrix *S,
								double *kl)
{
	let n = S0->size1;
	var L = gsl_matrix_alloc(n, n);
	var z = gsl_vector_alloc(n);
	var status = GSL_EDOM;

	if (!(observed_2ll > 0.0))
		goto done;

	gsl_matrix_memcpy(L, S0);
	if (gsl_linalg_cholesky_decomp1(L) != GSL_SUCCESS)
		goto done;

//	let dobjfn0 = sample_min2ll_from_cholesky(x->data, L);
//	printf("observed %f at start predicted dobjfn = %f\n", observed_2ll, dobjfn0);

	gsl_linalg_cholesky_solve(L, x, z);
	var a = 0.0;
	gsl_blas_ddot(x, z, &a);
	if (!(a > 0.0))
		goto done;

	/* Rank-one update S = S0 + delta*x*x'. By Sherman-Morrison,
	 * x'S^-1 x = a/(1 + delta*a), so choosing delta = 1/observed - 1/a
	 * makes it equal 'observed' exactly. S stays positive definite iff
	 * 1 + delta*a = a/observed > 0, which the checks above guarantee. */
	let delta = 1.0 / observed_2ll - 1.0 / a;
	forcount(i, n) {
		forcount(j, i + 1) {
			let xi = gsl_vector_get(x, i);
			let xj = gsl_vector_get(x, j);
			let v = gsl_matrix_get(S0, i, j) + delta * xi * xj;
			gsl_matrix_set(S, i, j, v);
			gsl_matrix_set(S, j, i, v);
		}
	}
	if (kl)
		*kl = 0.5 * (a / observed_2ll - 1.0 - log(a / observed_2ll));
	status = GSL_SUCCESS;

//	gsl_matrix_memcpy(L, S);
//	if (gsl_linalg_cholesky_decomp1(L) != GSL_SUCCESS)
//		goto done;
//	let dobjfn1 = sample_min2ll_from_cholesky(x->data, L);
//	printf("at end predicted dobjfn = %f\n", dobjfn1);

done:
	gsl_matrix_free(L);
	gsl_vector_free(z);
	return status;
}

/* Haar-random orthonormal basis in the columns of Q (n x n),
 * via modified Gram-Schmidt on Gaussian columns (applied twice for stability) */
static void random_orthonormal_basis(gsl_rng* rng, gsl_matrix* Q)
{
	let n = Q->size1;
	forcount(j, n) {
		var cj = gsl_matrix_column(Q, j);
		for (;;) {
			forcount(i, n)
				gsl_vector_set(&cj.vector, i, gsl_ran_ugaussian(rng));
			forcount(pass, 2) {
				forcount(k, j) {
					gsl_vector_view ck = gsl_matrix_column(Q, k);
					double d;
					gsl_blas_ddot(&ck.vector, &cj.vector, &d);
					gsl_blas_daxpy(-d, &ck.vector, &cj.vector);
				}
			}
			let nrm = gsl_blas_dnrm2(&cj.vector);
			if (nrm > 1e-8) {
				gsl_vector_scale(&cj.vector, 1. / nrm);
				break;
			}
		}
	}
}

/* observed dobjfn at a sample point; only logs to history when fitting */
static double resample_evaluate(const COVARGS* const covargs,
								const double* offset,
								const bool append_history)
{
	let objfn = encode_covariance_evaluate_offset(covargs->encode,
												  offset,
												  covargs->idata,
												  covargs->warm_etas,
												  covargs->advanfuncs,
												  append_history ? covargs->history : 0,
												  covargs->options);
	return objfn - covargs->f0;
}

static void resample_progress(const COVARGS* const covargs,
							  const int k,
							  const int nsample,
							  const double dobjfn0,
							  const double dobjfn1,
							  const double kl)
{
	char klmessage[256] = "";
	if (kl >= 0.)
		snprintf(klmessage, sizeof(klmessage), " kl %.6f", kl);
	char message[256];
	snprintf(message, sizeof(message), " sample %i/%i dobjfn (%.3f, %.3f)%s",
			 k+1, nsample, dobjfn0, dobjfn1, klmessage);
	var encode = covargs->encode;
	var extstream = covargs->extstream;
	var begin = covargs->begin;
	show_progress(&encode->popmodel, begin, covargs->details,
				  covargs->outstream, extstream, message);
}

#define RESAMPLE_MAX_INFLATE	10.	/* avg dobjfn floored at dobjfn0 / this */
#define RESAMPLE_CAP_ABOVE		6.	/* avg dobjfn capped at dobjfn0 + this */

static void print_resample_results(const double* predicted_dobjfn,
								   const double* observed_dobjfn,
								   const int nsample,
								   FILE* outstream,
								   FILE* diagstream)
{
	if (diagstream)
		fprintf(diagstream, OPENPMX_SFORMAT OPENPMX_SFORMAT OPENPMX_SFORMAT "\n",
				"SAMPLE", "PREDICTED", "OBSERVED");

	/* log ratio only defined for positive pairs; count the rest */
	var lratio = mallocvar(double, nsample);
	int nratio = 0;
	forcount(i, nsample) {
		let pred = predicted_dobjfn[i];
		let obs = observed_dobjfn[i];
		if (pred > 0. && obs > 0.)
			lratio[nratio++] = log(obs / pred);
		if (diagstream)
			fprintf(diagstream, OPENPMX_FFORMAT OPENPMX_FFORMAT OPENPMX_FFORMAT "\n",
					(double)i, pred, obs);
	}
	if (nratio > 1)
		info(outstream, "resample logratio %g(%g) negative %i\n",
			 gsl_stats_mean(lratio, 1, nratio),
			 gsl_stats_sd(lratio, 1, nratio),
			 nsample - nratio);
	else
		info(outstream, "resample logratio unavailable (%i negative)\n",
			 nsample - nratio);
	info(outstream, "resample predicted %g(%g) observed %g(%g)\n",
		 gsl_stats_mean(predicted_dobjfn, 1, nsample),
		 gsl_stats_sd(predicted_dobjfn, 1, nsample),
		 gsl_stats_mean(observed_dobjfn, 1, nsample),
		 gsl_stats_sd(observed_dobjfn, 1, nsample));
	free(lratio);
}

/* Additional samples for quadratic fitting: fixed Mahalanobis radius,
 * antithetic pairs, cov is updated after each pair and every evaluation
 * is logged to history. */
static void covariance_resample(const COVARGS* const covargs,
								gsl_matrix* cov,
								gsl_rng *rng,
								const int nsample_requested,
								const double expected_dobjfn)
{
	let nparam = (int)cov->size1;
	let outstream = covargs->outstream;

	/* whole blocks of nparam*2 so every forward sample has a backward partner */
	let blocksize = nparam * 2;
	let nblock = (nsample_requested + blocksize - 1) / blocksize;
	let nsample = nblock * blocksize;
	let radius = sqrt(expected_dobjfn);

	var predicted_dobjfn = mallocvar(double, nsample);
	var observed_dobjfn = mallocvar(double, nsample);
	var u = gsl_vector_alloc(nparam);
	var sample = gsl_vector_alloc(nparam);
	var L = gsl_matrix_alloc(nparam, nparam);
	var oldcov = gsl_matrix_alloc(nparam, nparam);
	var basis = gsl_matrix_alloc(nparam, nparam);
	gsl_matrix_set_identity(basis);	/* first block: identity */

	var k = 0;
	forcount(iblock, nblock) {
		/* L0 from cov at the START of the block, fixed for the whole block.
		 * Directions L0 q_j with Q orthonormal are conjugate, so the
		 * rank-one updates within the block do not interfere. */
		gsl_matrix_memcpy(L, cov);
		if (cholesky_decomposition(L)) {
			info(outstream, "error: resample cholesky failed\n");
			goto cleanup;
		}
		if (iblock > 0)
			random_orthonormal_basis(rng, basis);

		var block_kl = 0.;
		var block_capped = 0;
		var block_pairs = 0;
		var block_nonpos = 0;
		forcount(ibasis, nparam) {
			let q = gsl_matrix_column(basis, ibasis);

			/* forward sample, then its antithetic partner */
			var forward_dobjfn = 0.;
			forcount(backward, 2) {
				gsl_vector_memcpy(u, &q.vector);
				gsl_vector_scale(u, backward ? -radius : radius);
				gsl_blas_dgemv(CblasNoTrans, 1., L, u, 0., sample);

				let dobjfn0 = sample_min2ll_from_cholesky(sample->data, L);
				let dobjfn1 = resample_evaluate(covargs, sample->data, true);
				predicted_dobjfn[k] = dobjfn0;
				observed_dobjfn[k] = dobjfn1;

				var kl = -1.;
				if (!backward)
					forward_dobjfn = dobjfn1;	/* update once the pair is complete */
				else {
					let avg = (forward_dobjfn + dobjfn1) / 2.;
					if (!(avg > 0.)) {
						/* non-positive or NaN: a point better than f0, or a
						 * failed evaluation. No update possible. */
						block_nonpos += 1;
					} else {
						/* limit both shrinkage (cap) and inflation (floor) */
						let lo = dobjfn0 / RESAMPLE_MAX_INFLATE;
						let hi = dobjfn0 + RESAMPLE_CAP_ABOVE;
						let adj = fmin(fmax(avg, lo), hi);	/* avg == inf -> hi */
						gsl_matrix_memcpy(oldcov, cov);
						if (covariance_match_2ll(oldcov, sample, adj, cov, &kl) == GSL_SUCCESS) {
							block_pairs += 1;
							block_capped += (adj != avg);
							if (kl >= 0.)
								block_kl += kl;
						}
					}
				}
				resample_progress(covargs, k, nsample, dobjfn0, dobjfn1, kl);
				++k;
			}
		}
		info(outstream, "block %i/%i kl %.6f kl/nparam %.6f capped %i/%i nonpositive %i\n",
			 iblock+1, nblock, block_kl, block_kl / nparam,
			 block_capped, block_pairs, block_nonpos);
	}

	print_resample_results(predicted_dobjfn, observed_dobjfn, nsample, outstream, 0);

cleanup:
	gsl_matrix_free(basis);
	gsl_matrix_free(oldcov);
	gsl_matrix_free(L);
	gsl_vector_free(sample);
	gsl_vector_free(u);
	free(predicted_dobjfn);
	free(observed_dobjfn);
}

/* Diagnostics: chi(nparam) radius scaled to the target on average, cov is
 * not modified, nothing is logged to history. */
static void covariance_resample_diagnostics(const COVARGS* const covargs,
											FILE* diagstream,
											const gsl_matrix* const cov,
											gsl_rng *rng,
											const int nsample_requested,
											const double expected_dobjfn)
{
	let nparam = (int)cov->size1;
	let outstream = covargs->outstream;

	/* whole blocks of nparam*2 so every forward sample has a backward partner */
	let blocksize = nparam * 2;
	let nblock = (nsample_requested + blocksize - 1) / blocksize;
	let nsample = nblock * blocksize;

	/* cov is fixed, so one Cholesky factor serves every block */
	var L = gsl_matrix_alloc(nparam, nparam);
	gsl_matrix_memcpy(L, cov);
	if (cholesky_decomposition(L)) {
		info(outstream, "error: resample cholesky failed\n");
		gsl_matrix_free(L);
		return;
	}

	var predicted_dobjfn = mallocvar(double, nsample);
	var observed_dobjfn = mallocvar(double, nsample);
	var u = gsl_vector_alloc(nparam);
	var sample = gsl_vector_alloc(nparam);
	var basis = gsl_matrix_alloc(nparam, nparam);
	gsl_matrix_set_identity(basis);	/* first block: identity */

	var k = 0;
	forcount(iblock, nblock) {
		if (iblock > 0)
			random_orthonormal_basis(rng, basis);

		forcount(ibasis, nparam) {
			let q = gsl_matrix_column(basis, ibasis);
			let radius = sqrt(gsl_ran_chisq(rng, (double)nparam) * expected_dobjfn / nparam);

			/* forward sample, then its antithetic partner */
			forcount(backward, 2) {	
				gsl_vector_memcpy(u, &q.vector);
				gsl_vector_scale(u, backward ? -radius : radius);
				gsl_blas_dgemv(CblasNoTrans, 1., L, u, 0., sample);

				let dobjfn0 = sample_min2ll_from_cholesky(sample->data, L);
				let dobjfn1 = resample_evaluate(covargs, sample->data, false);
				predicted_dobjfn[k] = dobjfn0;
				observed_dobjfn[k] = dobjfn1;

				resample_progress(covargs, k, nsample, dobjfn0, dobjfn1, -1.);
				++k;
			}
		}
	}

	print_resample_results(predicted_dobjfn, observed_dobjfn, nsample, outstream, diagstream);

	gsl_matrix_free(basis);
	gsl_matrix_free(L);
	gsl_vector_free(sample);
	gsl_vector_free(u);
	free(predicted_dobjfn);
	free(observed_dobjfn);
}

typedef struct {
	double alpha;
	double z;
	double dobjfn;
} DELTA;

/* Critical p-value threshold
 * Split alpha in half to account for both the left and right tails
 * also calculate the delta objfn */
static DELTA dobjfn_from_alpha(const double alpha)
{
	let z = gsl_cdf_ugaussian_Qinv(alpha / 2.0);
	let dobjfn = z * z;
	return (DELTA) {
		.alpha = alpha,
		.z = z,
		.dobjfn = dobjfn,
	};
}

/* Fill in defaults for the user arguments. Returns the resolved
 * configuration and, if userargs is given, writes it back so the caller
 * can see what was actually used. */
static COVARIANCECONFIG covariance_args_resolve(COVARIANCECONFIG* const userargs, DELTA* delta)
{
	COVARIANCECONFIG args = { 0 };
	if (userargs) {
		args = *userargs;
		args.stage1 = stage1config_default(&userargs->stage1);
	}
	*delta = dobjfn_from_alpha(0.01);
	
	if (args.dobjfn <= 0.)
		args.dobjfn = delta->dobjfn;
	if (args.step_size <= 0.)
		args.step_size = 0.1;
	if (args.diagnostics.dobjfn <= 0.)
		args.diagnostics.dobjfn = args.dobjfn;
	if (userargs)
		*userargs = args;
	return args;
}

/* we have the upper and lower limits, so append to the .cov.ext file */
/* Value of the next limit entry if it matches (type, index), else 0.
 * Advances *k on a match. FIXED parameters have no entry and give 0. */
static double covext_limit_value(const COVLIMIT* const limit,
								 const int nlimit,
								 int* const k,
								 const int type,
								 const int index,
								 const bool upper)
{
	if (*k >= nlimit)
		return 0.;

	let ptr = &limit[*k];
	if ((int)ptr->param.type != type || ptr->param.index != index)
		return 0.;

	(*k)++;
	return upper ? ptr->upper : ptr->lower;
}

/* write one bound row (lower or upper) of the .cov.ext file */
static void covext_append_limit_row(FILE* extstream,
									const POPMODEL* const popmodel,
									const COVLIMIT* const limit,
									const int nlimit,
									const bool upper,
									const double runtime_s)
{
	fprintf(extstream, OPENPMX_IFORMAT,
			upper ? OPENPMX_COVFILE_UPPER_BOUND : OPENPMX_COVFILE_LOWER_BOUND);

	var k = 0;
	forcount(i, popmodel->ntheta) {
		let v = covext_limit_value(limit, nlimit, &k, PARAM_THETA, i, upper);
		fprintf(extstream, OPENPMX_TABLE_FORMAT, v);
	}
	forcount(i, popmodel->nsigma) {
		let v = covext_limit_value(limit, nlimit, &k, PARAM_SIGMA, i, upper);
		fprintf(extstream, OPENPMX_TABLE_FORMAT, v);
	}
	/* omega is written as the lower triangle, only diagonals have limits */
	forcount(i, popmodel->nomega) {
		forcount(j, i+1) {
			var v = 0.;
			if (i == j)
				v = covext_limit_value(limit, nlimit, &k, PARAM_OMEGA, i, upper);
			fprintf(extstream, OPENPMX_TABLE_FORMAT, v);
		}
	}
	assert(k == nlimit);

	fprintf(extstream, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(extstream, OPENPMX_TABLE_FORMAT, runtime_s);
	fprintf(extstream, OPENPMX_IFORMAT, 0);
	fprintf(extstream, "\n");
	fflush(extstream);
}

static COVLIMIT covlimit_make(ENCODE* const encode,
							  const gsl_matrix* const cov,
							  const int j,
							  const double z,
							  const int type,
							  const int index,
							  const double value,
							  const double* const decoded) /* points into encode->popmodel */
{
	decode_popmodel_limit(encode, cov, j, z);
	let upper = *decoded;
	decode_popmodel_limit(encode, cov, j, -z);
	let lower = *decoded;

	return (COVLIMIT) {
		.param = { .type = type,
				   .value = value,
				   .index = index },
		.lower = lower,
		.upper = upper,
	};
}

void pmx_covariance(OPENPMX* const source, COVARIANCECONFIG* const userargs)
{
	gsl_matrix *cov = 0;
	FILE* outstream = 0;
	FILE* extstream = 0;
	const ENCODELABEL* labels = 0;
	
	/* Work on a copy of the source object */
	var pmx = pmx_copy(source);
	pmx_ensure_state(&pmx);
	var pstate = pmx.state;

	/* open the results files */
	if (source->filename)
		outstream = estimate_covariance_results_fopen(source->filename, OPENPMX_COVFILE);
	info(outstream, "OpenPMX %i.%i.%i hash %s\n", OPENPMX_VERSION_MAJOR, OPENPMX_VERSION_MINOR, OPENPMX_VERSION_RELEASE, OPENPMX_GITHASH);
	if (source->filename)
		extstream = estimate_covariance_results_fopen(source->filename, OPENPMX_COVFILE OPENPMX_EXTFILE);
	if (!source->state) {
		info(outstream, "error: source does not have state, not previously estimated?\n");
		goto failed;
	}
	if (source->result.type == OBJFN_INVALID || source->result.objfn == DBL_MAX) {
		info(outstream, "error: source not previously estimated?\n");
		goto failed;
	}

	struct timespec begin;
	clock_gettime(CLOCK_MONOTONIC, &begin);

	/* copy the popmodel from the PMX object and make an encoder */
	ERRCTX errctx = { 0 };
	var popmodel = popmodel_init(pmx.theta, pmx.omega, pmx.sigma, &errctx);
	popmodel.result.objfn = source->result.objfn;
	if (errctx.len) 
		fatal(0, "%s: %s", __func__, errctx.errmsg);
	
	/* make an encoder for the popmodel */
	var encode = encode_init(&popmodel);
	encode_transform(&encode, &popmodel);
	let nparam = encode.nparam;
	if (nparam == 0) {
		info(outstream, "error: no parameters\n");
		goto failed;
	}
	labels = encode_labels_alloc(&encode, pmx.data._offset1);

	/* handle args and defaults */
	var options = options_init(&pmx);
	DELTA delta;
	let args = covariance_args_resolve(userargs, &delta);
	if (userargs)
		options.estimate.stage1 = args.stage1;
	if (!isfinite(args.dobjfn) || args.dobjfn <= 0.) {
		info(outstream, "error: invalid dobjfn (%g), covariance not calculated\n", args.dobjfn);
		goto failed;
	}

	info(outstream, "nparam %i initial step size %g\n", nparam, args.step_size);
	if (extstream)
		extfile_header(extstream, &popmodel, pmx.data._offset1, false);
	popmodel_eval_information(&popmodel, 0, 0, false, outstream, extstream, " final");

	/* do the evaluations */
	let idata = &pstate->idata;
	{
		/* allocate history if necessary */
		var history = (VB_HISTORY*)0;
		var _history = history_new(nparam);
		if (args.quadratic.nsample || args.type == COVARIANCE_SANDWICH_QUADRATIC) 
			history = &_history;

		/* center point f(0), needed by the adaptive step-size search */
		let f0 = source->result.objfn;
		var warm_etas = idata_etas_alloc(&source->state->idata);
		idata_etas_copy(&warm_etas, &source->state->idata);
		let covargs = (COVARGS) {
			.encode = &encode,
			.history = history,
			.idata = idata,
			.advanfuncs = pstate->advanfuncs,
			.options = &options,
			.details = args.details,
			.warm_etas = &warm_etas,
			.f0 = f0,
			.begin = &begin,
			.outstream = outstream,
			.extstream = extstream,
		};

		/* S matrix the OPG/BHHH cross-product-of-gradients */
		var S = (gsl_matrix*)0;
		if (args.type == COVARIANCE_S_MATRIX ||
			args.type == COVARIANCE_SANDWICH ||
			args.type == COVARIANCE_SANDWICH_QUADRATIC ) {
			S = pmx_opg_S(&covargs, labels, args.step_size);
		}
		
		/* Hessian matrix R */
		var R = (gsl_matrix*)0;
		if (args.type == COVARIANCE_R_MATRIX ||
			args.type == COVARIANCE_SANDWICH ||
			args.type == COVARIANCE_SANDWICH_QUADRATIC ) {
			R = pmx_hessian_R(&covargs, S, args.step_size);
		}

		/* examine S and R for deficiencies */
		if (S) 
			print_nonidentifiable_parameters(S, "S", labels, S_FLOOR_RELTOL, outstream);
		if (R) 
			print_nonidentifiable_parameters(R, "R", labels, R_FLOOR_RELTOL, outstream);

		/* make the covariance matrix via sandwich Rinv*S*Rinv, 2*Sinv, or 2*Rinv
		 * Or should it be 4*Sinv? */
		cov = gsl_matrix_alloc(nparam, nparam);
		if (!covariance_calculate(R, S, outstream, cov)) {
			gsl_matrix_free(cov);
			cov = 0;
			goto cov_cleanup;
		}

		/* make sure we have enough samples if we are doing quadratic smoothing */
		var nsample = args.quadratic.nsample;
		if (nsample < 0) {
			nsample = nparam * 2;
			let nfitparam = nparam * (nparam + 1) / 2;
			if (args.type == COVARIANCE_S_MATRIX) 
				nsample = nfitparam * 6;
			if (args.type == COVARIANCE_R_MATRIX) 
				nsample = nfitparam * 2;
		}

		/* additional resamples for quadratic smoothing */
		if (nsample > 0) {
			info(outstream, "additional samples %i\n", nsample);
			pmx_ensure_state_rng(source, &options);
			covariance_resample(&covargs,
								cov,
								source->state->rng,
								nsample,
								args.dobjfn);
		}
		
		/* quadratic smoothing */
		if (history) {
			info(outstream, "using quadratic smoothing\n");
			history_prune_above(history, f0 + 30., outstream);
			vb_fit_quadratic(history, f0, outstream, cov);
		}

		/* diagnostics if requested */
		if (args.diagnostics.nsample > 0) {
			FILE* diagstream = 0;
			if (args.diagnostics.filename)
				diagstream = fopen(args.diagnostics.filename, "w");

			let n = args.diagnostics.nsample;
			info(outstream, "diagnostics resample %i\n", n);
			pmx_ensure_state_rng(source, &options);
			covariance_resample_diagnostics(&covargs,
											diagstream,
											cov,
											source->state->rng,
											n,
											args.diagnostics.dobjfn);
			if (diagstream)
				fclose(diagstream);
		}

		/* cleanup */
cov_cleanup:
		idata_etas_free(&warm_etas);
		history_free(&_history);
		if (S)
			gsl_matrix_free(S);
		if (R)
			gsl_matrix_free(R);
	}
	
	/* did everything fail? */
	if (!cov) {
		info(outstream, "cov failed\n");
		goto failed;
	}

	/* This is dependent on the order of the parameter offsets in the
	 * ENCODE object. It would perhaps be better if there was a way to 
	 * indicate the relationship between the transformed and 
	 * untransformed spaces. */
	let z = delta.z;
	VECTOR(COVLIMIT) limit = { 0 };
	{
		var j = 0;
		vector_reserve(limit, nparam);
		forcount(i, popmodel.ntheta) {
			if (popmodel.thetaestim[i] == FIXED)
				continue;
			var l = covlimit_make(&encode,
								  cov,
								  j++,
								  z,
								  PARAM_THETA,
								  i,
								  popmodel.theta[i],
								  &encode.popmodel.theta[i]);
			vector_append(limit, l);
		}
		forcount(i, popmodel.nsigma) {
			if (popmodel.sigmafixed[i])
				continue;
			var l = covlimit_make(&encode,
								  cov,
								  j++,
								  z,
								  PARAM_SIGMA,
								  i,
								  popmodel.sigma[i],
								  &encode.popmodel.sigma[i]);
			vector_append(limit, l);
		}
		let rowcol = encode.omegainfo.nonfixed.rowcol;
		forcount(i, encode.omegainfo.nonfixed.n) {
			forcount(k, i+1) {
				let r = rowcol[i];
				let c = rowcol[k];
				if (popmodel.omegafixed[r][c])
					continue;
				/* off-diagonals get no limit (one encoded param can affect several) */
				if (i == k) {
					var l = covlimit_make(&encode,
										  cov,
										  j,
										  z,
										  PARAM_OMEGA,
										  r,
										  popmodel.omega[r][c],
										  &encode.popmodel.omega[r][c]);
					vector_append(limit, l);
				}
				j++;
			}
		}
		assert(j == encode.nparam);

		/* we have the upper and lower limits, so append to the .cov.ext file */
		struct timespec now;
		clock_gettime(CLOCK_MONOTONIC, &now);
		let runtime_s = timespec_time_difference(&begin, &now) / 1000.;

		covext_append_limit_row(extstream, &popmodel, limit.ptr, limit.size, false, runtime_s);
		covext_append_limit_row(extstream, &popmodel, limit.ptr, limit.size, true, runtime_s);
	}

	/* Eigendecomposition of the correlations in the approximate
	 * covariance matrix, the Hessian inverse. We normalize the diagonal
	 * to 1 to obtain the correlation matrix, so any scaling does not 
	 * play a role */
	var corr = gsl_matrix_alloc(nparam, nparam);
	{
		var Heigen = gsl_matrix_alloc(nparam, nparam);
		gsl_matrix_memcpy(Heigen, cov);
		scale_to_match_diagonal(Heigen, 0); /* scales for 1 on diagonal */
		gsl_matrix_memcpy(corr, Heigen);	/* save for printing correlations */

		var eval = gsl_vector_alloc(nparam);
		var evec = gsl_matrix_alloc(nparam, nparam);
		var w    = gsl_eigen_symmv_alloc(nparam);
		gsl_eigen_symmv(Heigen, eval, evec, w);
		gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);
		let lambda_min = gsl_vector_get(eval, 0);
		let lambda_max = gsl_vector_get(eval, nparam - 1);
		if (lambda_min != 0.)
			info(outstream, "condition number %g\n", fabs(lambda_max / lambda_min));
		else
			info(outstream, "cov is singular\n");

		gsl_eigen_symmv_free(w);
		gsl_matrix_free(evec);
		gsl_vector_free(eval);
		gsl_matrix_free(Heigen);
	}

	/* print requested upper and lower parameter limits */
	info(outstream, "dobjfn %g Z %g alpha %g\n", delta.dobjfn, delta.z, delta.alpha);
	info(outstream, OPENPMX_SFORMAT OPENPMX_SFORMAT OPENPMX_SFORMAT OPENPMX_SFORMAT "\n",
		"PARAMETER", "LOWER", "VALUE", "UPPER ");
	let offs = pmx.data._offset1 ? 1 : 0;
	forvector_ptr(p, limit) {
		char name[128];
		switch (p->param.type) {
		case PARAM_THETA:
			snprintf(name, sizeof(name), "THETA(%i)", p->param.index + offs);
			break;
		case PARAM_SIGMA:
			snprintf(name, sizeof(name), "SIGMA(%i)", p->param.index + offs);
			break;
		case PARAM_OMEGA:
			snprintf(name, sizeof(name), "OMEGA(%i,%i)", p->param.index + offs, p->param.index + offs);
			break;
		default:
			assert(0);
		}
		info(outstream, OPENPMX_SFORMAT OPENPMX_FFORMAT OPENPMX_FFORMAT OPENPMX_FFORMAT "\n", 
			name, p->lower, p->param.value, p->upper);
	}

	/* print correlations */
	print_correlations(corr, labels, outstream);
	gsl_matrix_free(corr);

	/* create the covariance object for the sources state */
	var covariance = &source->state->covariance;
	pmx_covariance_free(covariance);
	*covariance = (COVARIANCE) {
		.type = args.type,
		.limit = limit.ptr,
		.nlimit = limit.size,
		.dobjfn = args.dobjfn,
	};
	/* we DONT free the limit vector because we have taken the pointer
	 * into the covariance object. The memory gets freed when the 
	 * pmx_estimate() is repeated or the state is released */

	/* fallthrough to more cleanup */

failed:
    if (cov)
		gsl_matrix_free(cov);
    if (outstream)
		fclose(outstream);
	if (extstream)
		fclose(extstream);
	if (labels)
		encode_labels_free(labels);

	pmx_release_state(&pmx);
}

void pmx_covariance_free(COVARIANCE* cov)
{
	free((void*)cov->limit);
	cov->limit = 0;
	cov->nlimit = 0;
}

COVLIMIT pmx_covariance_limit(OPENPMX* const pmx, const PARAM* const param)
{
	let state = pmx->state;
	if (!state || state->covariance.limit == 0) 
		return (COVLIMIT) { 0 };

	const int offs = pmx->data._offset1 ? 1 : 0;
	let covariance = &state->covariance;
	forcount(i, covariance->nlimit) {
		let limit = &covariance->limit[i];
		let p = &limit->param;
		if (p->type == param->type && 
			p->index == param->index - offs) {

			/* internally, index is always 0-based, convert for user */
			var ret = covariance->limit[i];
			ret.param.index += offs;
			return ret;
		}
	}
	return (COVLIMIT) { 0 };
}
