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
 
/// This file calculates the covariance matrix of the estimated 
/// parameters. For the given alpha value (default 0.05) the upper and 
/// lower confidence limits are calculated. 

#include "openpmx.h"
#include "encode.h"
#include "print.h"
#include "scatter.h"
#include "stage1.h"
#include "pmxstate.h"
#include "defines.h"
#include "githash.h"
#include "linalg.h"
#include "utils/c22.h"
#include "utils/vector.h"
#include "utils/various.h"

#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_eigen.h>

/// This file calculates the covariance matrix of the estimate to provide
/// information about upper and lower bounds of paramater estimates.
///
/// Two matrices are computed, both by central finite differences: R,
/// the Hessian of the total objective function (NONMEM's R-matrix),
/// and S, the OPG/BHHH cross-product-of-gradients estimate (NONMEM's
/// S-matrix). These are combined into the robust sandwich estimator
/// Cov = Rinv * S * Rinv, which is what gets reported below.
///
/// The upper and lower bounds are not calculated for off-diagonal values
/// in omega. Its not really clear how to transform these from the
/// unbounded encoded space to the bounded paramater space.

/* 

Need to make a diagnostics script.
	For diagnostics see:
	https://chatgpt.com/share/6a8ab577-0768-83ed-bfee-af831ae642c9

 * */

#define HESSIAN_MSTEP_MAXITER				5
#define HESSIAN_FLOOR_RELTOL				1e-5

static void encode_covariance_evaluate(ENCODE* const test,
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

	popmodel->result.objfn = idata_objfn(idata, omegainfo->omega_nonzero_lndet);
	popmodel->result.type = OBJFN_CURRENT;
	popmodel->result.neval += 1;
	popmodel->result.nsig = 0.;
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
	var t = timespec_time_difference(begin, &now) / 1000.;
	popmodel_eval_information(popmodel, t, 0, details, outstream, extstream, extra);
}

typedef struct {
	ENCODE* const encode;
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

/* Perturbs along an arbitrary direction. For the S matrix the 
 * canonical-axis Jacobian is used but for the R matrix the steps are 
 * along the eigenvectors of S */
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
	var encode = covargs->encode;
	var idata = covargs->idata;
	var advanfuncs = covargs->advanfuncs;
	var options = covargs->options;
	var warm_etas = covargs->warm_etas;
	var f0 = covargs->f0;
	var begin = covargs->begin;
	var outstream = covargs->outstream;
	var extstream = covargs->extstream;
	let nparam = encode->nparam;
	let nindivid = idata->nindivid;

	var h_low = 0.;
	var h_high = 0.;
	var h_used = h;

	forcount(iter, HESSIAN_MSTEP_MAXITER) {
		OFFSET offset = { 0 };

		/* forward step */
		h_used = h;
		forcount(k, nparam) 
			offset[k] = h * direction[k];
		
		idata_etas_write(idata, warm_etas);
		encode_untransform(encode, offset);
		encode_covariance_evaluate(encode, idata, advanfuncs, options);
		char message2[256];
		snprintf(message2, sizeof(message2), "%s step %g", message, h_used);
		show_progress(&encode->popmodel, begin, covargs->details, outstream, extstream, message2);
		
		let fplus = encode->popmodel.result.objfn;
		if (iobjfn_plus_out)
			forcount(k, nindivid)
				iobjfn_plus_out[k] = idata->individ[k].iobjfn;

		/* backward step */
		forcount(k, nparam) 
			offset[k] = -h * direction[k];
			
		idata_etas_write(idata, warm_etas);
		encode_untransform(encode, offset);
		encode_covariance_evaluate(encode, idata, advanfuncs, options);
		show_progress(&encode->popmodel, begin, covargs->details, outstream, extstream, 0);
		
		let fminus = encode->popmodel.result.objfn;
		if (iobjfn_minus_out)
			forcount(k, nindivid)
				iobjfn_minus_out[k] = idata->individ[k].iobjfn;

		if (fplus_out)  
			*fplus_out  = fplus;
		if (fminus_out) 
			*fminus_out = fminus;
			
		/* check if average change falls within the target bracket */
		let upper = fabs(fplus - f0);
		let lower = fabs(fminus - f0);
		let avgdiff = (upper + lower) / 2.;
		if (avgdiff > target_delta_lower && 
			avgdiff < target_delta_upper) {
			return h_used;
		}

		/* dobjfn scales as ~0.5*kappa*h^2 locally, so rather than
		 * bisecting, solve directly for the h expected to hit the
		 * target. Converges in far fewer evaluations than doubling/
		 * bisection because it uses the known h^2 scaling instead
		 * of just monotonicity. */
		let target = sqrt(target_delta_lower * target_delta_upper);
		var ratio = (avgdiff > 0.) ? sqrt(target / avgdiff) : 4.;
		ratio = fmin(fmax(ratio, 0.1), 10.); /* guard a wild jump off one noisy sample */
		var hnext = h * ratio;

		/* h_low/h_high kept only as a fallback bracket for cases where
		 * the local quadratic model fails (e.g. BLQ discontinuities
		 * making dobjfn(h) non-monotonic) */
		if (avgdiff <= target_delta_lower)
			h_low = h;
		else
			h_high = h;
		if (h_low > 0. && h_high > 0. && (hnext <= h_low || hnext >= h_high))
			hnext = (h_low + h_high) / 2.; /* jump left the bracket, fall back to bisection */

		h = hnext;
	}
		
	info(outstream, "cov: warning: eigenmode step size did not reach objfn "
		"change of [%g,%g] after %i iterations, using h=%g\n",
		target_delta_lower, target_delta_upper, 
		HESSIAN_MSTEP_MAXITER, h_used);
		
	return h_used;
}

/* R matrix: the Hessian of the *total* objective function (summed
 * over all individuals), via central second differences in the
 * transformed parameter space. This is distinct from the OPG/BHHH
 * "S" matrix (S/cov below, built from per-individual score outer
 * products) -- R uses second derivatives of the objfn directly, no
 * per-individual decomposition needed. Combining both gives the
 * robust sandwich estimator Cov = Rinv * S * Rinv. */
static gsl_matrix* pmx_hessian_R(const COVARGS* const covargs,
								 const gsl_matrix* const S)
{
	var encode = covargs->encode;
	var idata = covargs->idata;
	var advanfuncs = covargs->advanfuncs;
	var options = covargs->options;
	var warm_etas = covargs->warm_etas;
	var f0 = covargs->f0;
	var begin = covargs->begin;
	var outstream = covargs->outstream;
	var extstream = covargs->extstream;

	let nparam = encode->nparam;

	/* rotate: eigendecompose S to get the directions to step along.
	 * Scopy protects the caller's S from gsl_eigen_symmv's in-place
	 * destruction -- S is still needed untouched for the sandwich
	 * estimator after this returns. */
	var evalS = gsl_vector_alloc(nparam);
	var U     = gsl_matrix_alloc(nparam, nparam);
	{
		var Scopy = gsl_matrix_alloc(nparam, nparam);
		gsl_matrix_memcpy(Scopy, S);
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
	let COVARIANCE_DOBJFN_TARGET_LOWER = 2;
	let COVARIANCE_DOBJFN_TARGET_UPPER = 4.;
	let avg_delta = (COVARIANCE_DOBJFN_TARGET_LOWER + COVARIANCE_DOBJFN_TARGET_UPPER) / 2.;
	typeof(encode->offset) step_size; /* per eigenmode, filled in when doing diagonal */
	{
		forcount(p, nparam) {
			OFFSET fplus;
			OFFSET fminus;
			OFFSET direction;
			forcount(k, nparam) 
				direction[k] = gsl_matrix_get(U, k, p);

			char message[128];
			snprintf(message, sizeof(message), " R(%i,%i)", p, p);

			let lam = fmax(gsl_vector_get(evalS, p), HESSIAN_FLOOR_RELTOL * lambda_max);
			let h_init = 2. * sqrt(avg_delta) / sqrt(lam);
			step_size[p] = pmx_hessian_stepsize_dir(covargs, direction, h_init,
													0, 0, 
													&fplus[p], &fminus[p],
													COVARIANCE_DOBJFN_TARGET_LOWER, 
													COVARIANCE_DOBJFN_TARGET_UPPER,
													message);
			let h = step_size[p];
			let hpp = (fplus[p] - 2. * f0 + fminus[p]) / (h * h);
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
		
			/* both positive */
			OFFSET offset = { 0 };
			forcount(k, nparam)
				offset[k] = hp * gsl_matrix_get(U, k, p) + hq * gsl_matrix_get(U, k, q);

			idata_etas_write(idata, warm_etas);
			encode_untransform(encode, offset);
			encode_covariance_evaluate(encode, idata, advanfuncs, options);
			let fpp = encode->popmodel.result.objfn;
			char message[128];
			snprintf(message, sizeof(message), " R(%i,%i)", p, q);
			show_progress(&encode->popmodel, begin, covargs->details, outstream, extstream, message);

			/* one positive, one negative */
			forcount(k, nparam)
				offset[k] = hp * gsl_matrix_get(U, k, p) - hq * gsl_matrix_get(U, k, q);

			idata_etas_write(idata, warm_etas);
			encode_untransform(encode, offset);
			encode_covariance_evaluate(encode, idata, advanfuncs, options);
			let fpm = encode->popmodel.result.objfn;
			show_progress(&encode->popmodel, begin, covargs->details, outstream, extstream, 0);

			/* one negative, one positive */
			forcount(k, nparam)
				offset[k] = -hp * gsl_matrix_get(U, k, p) + hq * gsl_matrix_get(U, k, q);

			idata_etas_write(idata, warm_etas);
			encode_untransform(encode, offset);
			encode_covariance_evaluate(encode, idata, advanfuncs, options);
			let fmp = encode->popmodel.result.objfn;
			show_progress(&encode->popmodel, begin, covargs->details, outstream, extstream, 0);

			/* both negative */
			forcount(k, nparam)
				offset[k] = -hp * gsl_matrix_get(U, k, p) - hq * gsl_matrix_get(U, k, q);

			idata_etas_write(idata, warm_etas);
			encode_untransform(encode, offset);
			encode_covariance_evaluate(encode, idata, advanfuncs, options);
			let fmm = encode->popmodel.result.objfn;
			show_progress(&encode->popmodel, begin, covargs->details, outstream, extstream, 0);

			/* set off-diagonal of symmetric R matrix based on samples */
			let hpq = (fpp - fpm - fmp + fmm) / (4. * hp * hq);
			gsl_matrix_set(Rp, p, q, hpq);
			gsl_matrix_set(Rp, q, p, hpq);
		}
	}

	/* reset the old eta state */
	idata_etas_write(idata, warm_etas);

	/* de-rotate: R = U * Rp * U' back to canonical parameter basis --
	 * from here on R is an ordinary Hessian in the same basis as S,
	 * usable by the existing regularize/sandwich/report code unchanged */
	var R = gsl_matrix_alloc(nparam, nparam);
	{
		var tmp = gsl_matrix_alloc(nparam, nparam);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., U, Rp, 0., tmp);
		gsl_blas_dgemm(CblasNoTrans, CblasTrans,   1., tmp, U, 0., R);
		gsl_matrix_free(tmp);
	}
	
	gsl_matrix_free(Rp);
	gsl_matrix_free(U);

	return R;
}

/* Invert a symmetric matrix via eigendecomposition, flooring any
 * eigenvalue below reltol * lambda_max before inverting. Finite-
 * difference noise can push R's smallest curvature directions to
 * zero or slightly negative; this regularizes those away instead of
 * requiring R to be exactly positive definite for a plain Cholesky
 * inversion to succeed. Returns the number of eigenvalues that were
 * floored (0 = ordinary PD inverse, identical to what Cholesky would
 * give), or -1 if R has no positive curvature at all (lambda_max <= 0),
 * which is a genuine degeneracy this cannot regularize away. */
static int spd_regularize_invert(gsl_matrix* const Minv,
								 const gsl_matrix* const M,
								 const gsl_matrix* const P,
								 const double reltol)
{
	let n = (int)M->size1;

	var A = gsl_matrix_alloc(n, n);
	gsl_matrix_memcpy(A, M);
	forcount(i, n) {
		for (int k = i + 1; k < n; ++k) {
			let avg = 0.5 * (gsl_matrix_get(A, i, k) + gsl_matrix_get(A, k, i));
			gsl_matrix_set(A, i, k, avg);
			gsl_matrix_set(A, k, i, avg);
		}
	}

	var eval = gsl_vector_alloc(n);
	var evec = gsl_matrix_alloc(n, n);
	var w = gsl_eigen_symmv_alloc(n);
	gsl_eigen_symmv(A, eval, evec, w);
	gsl_eigen_symmv_free(w);
	gsl_matrix_free(A);

	let lambda_max = gsl_vector_max(eval);
	if (lambda_max <= 0.) {
		gsl_matrix_free(evec);
		gsl_vector_free(eval);
		return -1;
	}

	/* Stability — use S as an informed regularization prior, not an
	 * arbitrary floor */
	/* https://claude.ai/share/8a8cca70-6c13-4efb-9ef4-57487e02b719 */
	let floor = reltol * lambda_max;
	var Pv = gsl_vector_alloc(n);
	int nfloored = 0;
	forcount(i, n) {
		if (gsl_vector_get(eval, i) < floor) {
			var vi = gsl_matrix_column(evec, i);
			gsl_blas_dgemv(CblasNoTrans, 1., P, &vi.vector, 0., Pv);
			double lambda_prior;
			gsl_blas_ddot(&vi.vector, Pv, &lambda_prior);
			gsl_vector_set(eval, i, fmax(0.5 * lambda_prior, floor));
			++nfloored;
		}
	}
	gsl_vector_free(Pv);
    
	/* Minv = V * diag(1/eval) * V' */
	var Vscaled = gsl_matrix_alloc(n, n);
	gsl_matrix_memcpy(Vscaled, evec);
	forcount(k, n) {
		var col = gsl_matrix_column(Vscaled, k);
		gsl_vector_scale(&col.vector, 1. / gsl_vector_get(eval, k));
	}
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., Vscaled, evec, 0., Minv);

	gsl_matrix_free(Vscaled);
	gsl_matrix_free(evec);
	gsl_vector_free(eval);

	return nfloored;
}



/* S matrix: the OPG/BHHH cross-product-of-gradients estimate,
 * S = J'J, where J is the per-individual Jacobian of iobjfn with
 * respect to each transformed parameter (central differences).
 * S is kept uninverted here, same as pmx_hessian_R above, so the caller
 * decides whether to invert alone or combine into a sandwich
 * estimator with R. */
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

	/* Single Pass: Compute OPG Jacobian along canonical coordinate axes
	 * also known as the BHHH or OPG/Fisher information estimate */
	let COVARIANCE_DOBJFN_TARGET_LOWER = 0.5;
	let COVARIANCE_DOBJFN_TARGET_UPPER = 2.;
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

	typeof(encode->offset) offset = { 0 };
	offset[j] = z * SD;
	
	encode_untransform(encode, offset);
}

static gsl_matrix* covariance_alloc(const gsl_matrix* const R, 
									const gsl_matrix* const S,
									FILE* outstream)
{
	let nparam = S->size2;
	var cov = gsl_matrix_alloc(nparam, nparam);

	/* Sandwich estimator: cov = Rinv * S * Rinv */
	if (R) {
		info(outstream, "cov using Sandwich\n");

		var Rinv = gsl_matrix_alloc(nparam, nparam);
		gsl_matrix_memcpy(Rinv, R);
		if (gsl_linalg_cholesky_decomp(Rinv) != GSL_SUCCESS ||
			gsl_linalg_cholesky_invert(Rinv) != GSL_SUCCESS) {
			info(outstream, "cov: warning: singular or ill-conditioned Hessian (R) matrix\n");

			/* cholesky failed, do a safer inversion */
			let nfloored = spd_regularize_invert(Rinv, R, S, HESSIAN_FLOOR_RELTOL);
			if (nfloored < 0) {
				info(outstream, "cov: error: Hessian (R) has no positive curvature in any direction\n");
				gsl_matrix_free(Rinv);
				goto cov_from_S_matrix;
			}
			if (nfloored > 0) {
				info(outstream, "cov: warning: %i/%i eigenvalue(s) of R were <= %g x largest "
					"(likely finite-difference noise), floored before inversion\n",
					nfloored, nparam, HESSIAN_FLOOR_RELTOL);
			}
		}

		var tmp = gsl_matrix_alloc(nparam, nparam);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Rinv, S, 0.0, tmp);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmp, Rinv, 0.0, cov);
		gsl_matrix_free(tmp);
		gsl_matrix_free(Rinv);
	
	/* S estimator: Cov = 2*Sinv */
	} else {

cov_from_S_matrix:
		info(outstream, "cov using S matrix\n");

		gsl_matrix_memcpy(cov, S);
		if (gsl_linalg_cholesky_decomp(cov) != GSL_SUCCESS ||
			gsl_linalg_cholesky_invert(cov) != GSL_SUCCESS) {
			info(outstream, "cov: error: singular or ill-conditioned OPG/BHHH matrix\n");
			gsl_matrix_free(cov);
			return 0;
		}
		gsl_matrix_scale(cov, 2.);
	}
	return cov;
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
	VECTOR(CORRENTRY) entries = { 0 };
	
	let nparam = corr->size1;
	vector_reserve(entries, nparam * (nparam + 1) / 2);
	forcount(i, nparam) {
		forcount(j, i) {
			vector_append(entries, 
				(CORRENTRY) {
					.i = labels[i],
					.j = labels[j],
					.val = gsl_matrix_get(corr, i, j),
				}
			);
		}
	}
	qsort(entries.mutptr, entries.size, sizeof(entries.mutptr[0]), correntry_sort_callback);

	info(outstream, "cov correlations\n");
	info(outstream, OPENPMX_SFORMAT OPENPMX_SFORMAT OPENPMX_SFORMAT "\n", "CORRELATION", "PARAMETER1", "PARAMETER2");
	forvector_ptr(p, entries) 
		info(outstream, OPENPMX_FFORMAT OPENPMX_SFORMAT OPENPMX_SFORMAT "\n", p->val, p->i, p->j);

	vector_free(entries);
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

/* Finds ALL non-identifiable parameters from a raw (unnormalized) S
 * matrix -- not just the single smallest eigenvalue. Sums each
 * parameter's squared component across every eigenvalue <= reltol *
 * lambda_max, i.e. the diagonal of the projector onto that whole
 * near-null subspace. This is basis-independent: within a degenerate
 * or near-degenerate cluster of small eigenvalues, individual
 * eigenvectors are only defined up to rotation, but the projector
 * onto the subspace they jointly span is not -- so this correctly
 * finds multiple simultaneous non-identifiabilities (a fully unused
 * parameter, two parameters only identifiable as a sum, etc.)
 * without depending on which particular basis GSL happened to return.
 * Score is in [0,1]; near 1 means "this parameter's information is
 * (almost) entirely confined to the rank-deficient subspace". */
static void print_nonidentifiable_parameters(const gsl_matrix* const S,
											 const ENCODELABEL* labels,
											 const double reltol,
											 FILE* outstream)
{
	let nparam = (int)S->size1;

	var Scopy = gsl_matrix_alloc(nparam, nparam);
	gsl_matrix_memcpy(Scopy, S);
	var eval = gsl_vector_alloc(nparam);
	var evec = gsl_matrix_alloc(nparam, nparam);
	var w    = gsl_eigen_symmv_alloc(nparam);
	gsl_eigen_symmv(Scopy, eval, evec, w);
	gsl_eigen_symmv_free(w);
	gsl_matrix_free(Scopy);

	let lambda_max = gsl_vector_max(eval);
	let floor = reltol * lambda_max;

	VECTOR(EIGCONTRIB) entries = { 0 };
	vector_reserve(entries, nparam);
	forcount(i, nparam) 
		vector_append(entries, (EIGCONTRIB) { .name = labels[i], .contrib = 0. });

	int nnull = 0;
	forcount(k, nparam) {
		if (gsl_vector_get(eval, k) > floor)
			continue;
		++nnull;
		forcount(i, nparam) {
			let v = gsl_matrix_get(evec, i, k);
			entries.mutptr[i].contrib += v * v;
		}
	}
	gsl_matrix_free(evec);
	gsl_vector_free(eval);

	qsort(entries.mutptr, entries.size, sizeof(entries.mutptr[0]), eigcontrib_sort_callback);

	info(outstream, "cov identifiability %i of %i eigenvalue(s) <= %g x largest\n",
		nnull, nparam, reltol);
	if (nnull) {
		info(outstream, OPENPMX_SFORMAT OPENPMX_SFORMAT "\n", "PARAMETER", "NULLSPACE_CONTRIB");
		forvector_ptr(p, entries)
			info(outstream, OPENPMX_SFORMAT OPENPMX_FFORMAT "\n", p->name, p->contrib);
	}
	
	vector_free(entries);
}

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

static void diagnostics_nsample(const COVARGS* const covargs,
								const gsl_matrix* cov,
								gsl_rng *rng,
								const int nsample)
{
	let nparam = (int)cov->size1;
	let outstream = covargs->outstream;

	info(outstream, "cov nsample %i\n", nsample);
	
	var L = gsl_matrix_alloc(nparam, nparam);
	gsl_matrix_memcpy(L, cov);
	let err = cholesky_decomposition(L);
	if (err) {
		info(outstream, "error: cov: nsample: Cholesky of covariance failed\n");
		goto failed;
	}

	var predicted_dobjfn = mallocvar(double, nsample);
	var observed_dobjfn = mallocvar(double, nsample);

    var u = gsl_vector_alloc(nparam);
    var sample = gsl_vector_alloc(nparam);
	forcount(k, nsample) {
		var encode = covargs->encode;
		var idata = covargs->idata;
		var advanfuncs = covargs->advanfuncs;
		var options = covargs->options;
		var warm_etas = covargs->warm_etas;
		var begin = covargs->begin;
		var extstream = covargs->extstream;
		var f0 = covargs->f0;

		/* diagonal sample */
		forcount(i, nparam)
			gsl_vector_set(u, i, gsl_ran_ugaussian(rng));

		/* rotate to covariance using Cholesky */
		gsl_blas_dgemv(CblasNoTrans, 1., L, u, 0., sample);

		/* expected objfn */
		let dobjfn0 = sample_min2ll_from_cholesky(sample->data, L);

		/* actually do the sample */
		idata_etas_write(idata, warm_etas);
		encode_untransform(encode, sample->data);
		encode_covariance_evaluate(encode, idata, advanfuncs, options);
		let dobjfn1 = encode->popmodel.result.objfn - f0;
		char message[256];
		snprintf(message, sizeof(message),
			" sample %i/%i dobjfn (%.3f, %.3f)",
			k, nsample, dobjfn0, dobjfn1);
		show_progress(&encode->popmodel, begin, covargs->details, outstream, extstream, message);

		/* save results */
		predicted_dobjfn[k] = dobjfn0;
		observed_dobjfn[k] = dobjfn1;
	}

	/* print results */
	info(outstream, "cov nsample\n");
	var ratio = mallocvar(double, nsample);
	info(outstream, OPENPMX_SFORMAT OPENPMX_SFORMAT OPENPMX_SFORMAT OPENPMX_SFORMAT "\n",
		"SAMPLE", "PREDICTED", "OBSERVED", "RATIO");
	forcount(k, nsample) {
		let pred = predicted_dobjfn[k];
		let obs = observed_dobjfn[k];
		ratio[k] = obs / pred;
		info(outstream, OPENPMX_IFORMAT OPENPMX_FFORMAT OPENPMX_FFORMAT OPENPMX_FFORMAT "\n",
			k, pred, obs, ratio[k]);
	}
	info(outstream, "cov nsample ratio mean %g sd %g\n",
		 gsl_stats_mean(ratio, 1, nsample),
		 gsl_stats_sd(ratio, 1, nsample));

	gsl_vector_free(sample);
	gsl_vector_free(u);
	free(ratio);
	free(predicted_dobjfn);
	free(observed_dobjfn);

failed:
	gsl_matrix_free(L);
}

void pmx_covariance(OPENPMX* const source, COVARIANCECONFIG* const userargs)
{
	gsl_matrix *S = 0;
	gsl_matrix *R = 0;
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
		info(outstream, "cov: error: source does not have state, not previously estimated?\n");
		goto failed;
	}

	/* handle default args */
	var options = options_init(&pmx);
	COVARIANCECONFIG args = { 0 };
	if (userargs) {
		args = *userargs;
		args.stage1 = stage1config_default(&userargs->stage1);
		options.estimate.stage1 = args.stage1; 
	}
	if (args.alpha <= 0.)
		args.alpha = 0.05;
	if (args.step_size <= 0.)
		args.step_size = options.estimate.step_refine;
	if (userargs)
		*userargs = args;

	struct timespec begin;
	clock_gettime(CLOCK_MONOTONIC, &begin);

	/* copy the popmodel from the PMX object and make an encoder */
	ERRCTX errctx = { 0 };
	let popmodel = popmodel_init(pmx.theta, pmx.omega, pmx.sigma, &errctx);
	if (errctx.len) 
		fatal(0, "cov: %s: %s", __func__, errctx.errmsg);

	var encode = encode_init(&popmodel);
	encode_transform(&encode, &popmodel);
	let nparam = encode.nparam;
	labels = encode_labels_alloc(&encode, pmx.data._offset1);

	info(outstream, "cov nparam %i step_size %g\n", nparam, args.step_size);
	if (extstream)
		extfile_header(extstream, &encode.popmodel, pmx.data._offset1);

	if (!isfinite(args.alpha) || args.alpha <= 0. || args.alpha >= 1.) {
		info(outstream, "cov: error: invalid alpha (%g), covariance not calculated\n", args.alpha);
		goto failed;
	}
	if (nparam == 0) {
		info(outstream, "cov: error: no estimated parameters, covariance not calculated\n");
		goto failed;
	}
	let idata = &pstate->idata;
	let nindivid = idata->nindivid;
	if (nindivid < nparam) {
		info(outstream, "cov: error: nindivid < nparam, covariance not calculated\n");
		goto failed;
	}

	/* do the evaluations */
	{
		/* center point f(0), needed by the adaptive step-size search */
		var warm_etas = idata_etas_alloc(&source->state->idata);
		idata_etas_copy(&warm_etas, &source->state->idata);
		let f0 = source->result.objfn;
		let covargs = (COVARGS) {
			.encode = &encode,
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

		/* OPG/BHHH cross-product-of-gradients matrix S */
		info(outstream, "cov S\n");
		S = pmx_opg_S(&covargs, labels, args.step_size);

		/* diagnostics for the S matrix */
		print_nonidentifiable_parameters(S, labels, HESSIAN_FLOOR_RELTOL, outstream);

		/* Hessian matrix R */
		if (args.type == COVARIANCE_SANDWICH) {
			info(outstream, "cov R\n");
			R = pmx_hessian_R(&covargs, S);
		}
		
		/* make the covariance matrix via Sandwich RinvT*S*Rinv or 2*Sinv */
		cov = covariance_alloc(R, S, outstream);

		/* resample if requested */
		if (args.nsample) {
			pmx_ensure_state_rng(source, &options);
			diagnostics_nsample(&covargs, cov, source->state->rng, args.nsample);
		}

		/* cleanup */
		idata_etas_free(&warm_etas);
	}
	
	/* did everything fail? */
	if (!cov) {
		info(outstream, "cov failed\n");
		goto failed;
	}

	/* Critical p-value threshold
	 * Split alpha in half to account for both the left and right tails */
	let z = gsl_cdf_ugaussian_Qinv(args.alpha / 2.0);

	/* This is dependent on the order of the parameter offsets in the
	 * ENCODE object. It would perhaps be better if there was a way to 
	 * indicate the relationship between the transformed and 
	 * untransformed spaces */
	VECTOR(COVLIMIT) limit = { 0 };
	var corr = gsl_matrix_alloc(nparam, nparam);
	{
		var j = 0;
		vector_reserve(limit, nparam);
		forcount(i, popmodel.ntheta) {
			if (popmodel.thetaestim[i] != FIXED) {
				decode_popmodel_limit(&encode, cov, j, z);
				let upperv = encode.popmodel.theta[i];

				decode_popmodel_limit(&encode, cov, j, -z); 
				let lowerv = encode.popmodel.theta[i];

				/* construct the COVLIMIT and append */
				var covlimit = (COVLIMIT) {
					.param = {
						.type = PARAM_THETA,
						.value = popmodel.theta[i],
					},
					.lower = lowerv,
					.upper = upperv,
				};
				vector_append(limit, covlimit);

				++j;
			}
		}
		
		forcount(i, popmodel.nsigma) {
			if (popmodel.sigmafixed[i] == 0) {
				decode_popmodel_limit(&encode, cov, j, z);
				let upperv = encode.popmodel.sigma[i];

				decode_popmodel_limit(&encode, cov, j, -z); 
				let lowerv = encode.popmodel.sigma[i];

				/* construct the COVLIMIT and append */
				var covlimit = (COVLIMIT) {
					.param = {
						.type = PARAM_SIGMA,
						.value = popmodel.sigma[i],
					},
					.lower = lowerv,
					.upper = upperv,
				};
				vector_append(limit, covlimit);

				++j;
			}
		}

		forcount(i,  encode.omegainfo.nonfixed.n) {
			forcount(k, i+1) {
				var rowcol = encode.omegainfo.nonfixed.rowcol;
				let r = rowcol[i];
				let c = rowcol[k];
				if (popmodel.omegafixed[r][c] == 0) {
					if (i == k) {
						decode_popmodel_limit(&encode, cov, j, z);
						let upperv = encode.popmodel.omega[r][c];

						decode_popmodel_limit(&encode, cov, j, -z); 
						let lowerv = encode.popmodel.omega[r][c];

						/* construct the COVLIMIT and append */
						var covlimit = (COVLIMIT) {
							.param = {
								.type = PARAM_OMEGA,
								.value = popmodel.omega[r][c],
							},
							.lower = lowerv,
							.upper = upperv,
						};
						vector_append(limit, covlimit);

					} else {
						/* since encoding of off diagonals could probably
						 * influence more than one estimated parameter then
						 * its difficult to know the best way of reporting
						 * this. So for now we simply dont calculate SD
						 * for omega off-diagonals */
					}

					++j;
				}
			}
		}
		assert(j == encode.nparam);
	}

	/* Eigendecomposition of the correlations in the approximate
	 * covariance matrix, the Hessian inverse. We normalize the diagonal
	 * to 1 to obtain the correlation matrix, so any scaling does not 
	 * play a role */
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
			info(outstream, "cov condition number %g\n", fabs(lambda_max / lambda_min));
		else
			info(outstream, "cov is singular\n");

		gsl_eigen_symmv_free(w);
		gsl_matrix_free(evec);
		gsl_vector_free(eval);
		gsl_matrix_free(Heigen);
	}

	/* print requested upper and lower parameter limits */
	info(outstream, "cov alpha %g Z %g\n", args.alpha, z);
	info(outstream, OPENPMX_SFORMAT OPENPMX_SFORMAT OPENPMX_SFORMAT OPENPMX_SFORMAT "\n",
		"PARAMETER", "LOWER", "VALUE", "UPPER ");
	forvector_ptr(p, limit) {
		char name[128];
		switch (p->param.type) {
		case PARAM_THETA:
			snprintf(name, sizeof(name), "THETA(%i)", p->param.index);
			break;
		case PARAM_SIGMA:
			snprintf(name, sizeof(name), "SIGMA(%i)", p->param.index);
			break;
		case PARAM_OMEGA:
			snprintf(name, sizeof(name), "OMEGA(%i)", p->param.index);
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
	covariance_free(covariance);
	*covariance = (COVARIANCE) {
		.type = args.type,
		.limit = limit.ptr,
		.nlimit = limit.size,
		.alpha = args.alpha,
	};
	/* we DONT free the limit vector because we have taken the pointer
	 * into the covariance object. The memory gets freed when the 
	 * pmx_estimate() is repeated or the state is released */

	/* fallthrough to more cleanup */

failed:
	if (S)
		gsl_matrix_free(S);
	if (R)
		gsl_matrix_free(R);
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

void covariance_free(COVARIANCE* cov)
{
	free((void*)cov->limit);
}

COVLIMIT pmx_covariance_limit(OPENPMX* const pmx, const PARAMETER* const param)
{
	let state = pmx->state;
	if (!state || state->covariance.limit == 0) 
		fatal(0, "cov: error: source does not have covariance limits\n");

	let covariance = &state->covariance;
	forcount(i, covariance->nlimit) {
		let limit = &covariance->limit[i];
		let p = &limit->param;
		if (p->type == param->type && 
			p->index == param->index) 
			return covariance->limit[i];
	}
	return (COVLIMIT) { 0 };
}

