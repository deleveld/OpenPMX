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

/// This advancer uses eigendecomposition to obtain exact analytical
/// (no step size error) solutions for any linear compartmental model,
/// regardless of the number of compartments or topology.
///
/// Instead of hardcoding the analytical solution for a specific model
/// structure or using numerical ODE integration, this advancer:
///
///   1. User must supply system matrix
///   2. Eigendecomposes the matrix: A = V * diag(eigvals) * V^-1
///   3. Advances state analytically in modal coordinates
///
/// This is equivalent to NONMEM's ADVAN5 (linear general compartmental
/// model). The technique is described in:
///
/// - Joachim J, "TCI Eigendecomposition" (2026)
/// - Minto C, Schnider T, "From Eigenvalues to Plasma Coefficients" (2026)
///
/// The modal state update for constant input over interval dt:
///
///   x_modal(t+dt) = x_modal(t) * exp(eigval*dt)
///                  + (V^-1 * J) / (-eigval) * (1 - exp(eigval*dt))
///
/// where eigval < 0 for stable systems, V is the eigenvector matrix,
/// and J is the input vector (infusion rates into each compartment).
///
/// Uses libGSL for eigendecomposition (gsl_eigen_nonsymmv) and matrix
/// inversion (LU decomposition).

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "advan/advan.h"
#include "utils/c22.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

/// For the eigensystem advancer the system matrix must be provided
/// in row-major order. The user can write this using
/// `pmx_advan_eigen_sysmat()` or via openpmxtran as `SYSMAT()`.
/*------------------------------------------------------------------------
 * sysmat is stored in row-major order (C convention, GSL convention):
 *   sysmat[i * n + j] = element at row i, column j
 *
 * The matrix operates on amounts (consistent with OpenPMX state).
 * Infusion rates are handled separately by the advancer.
 *----------------------------------------------------------------------*/

typedef struct {
	ADVAN advan;

	/* system matrix */
	/* user will fill this in in $IMODEL() call, via SYSMAT() macro */
	/* cached sysmat to avoid unnecessary recalculation */
	double* sysmat_data;
	double* last_sysmat_data;
	int last_recalc_initcount;
	
	/* Cached eigendecomposition results */
	int n;										/* matrix dimension */
	
	/* GSL workspaces allocated once in constructor, reused */
	gsl_eigen_nonsymmv_workspace* eigen_ws;
	gsl_vector_complex* eval;
	gsl_matrix_complex* evec;
	gsl_permutation* perm;    
	double* eigvals;							/* eigenvalues (real parts) */
	
	gsl_matrix* V_m;							/* right eigenvector matrix (tda == n, contiguous row-major) */
	gsl_matrix* Vinv_m;							/* inverse of V (tda == n, contiguous row-major) */

	/* transposed copies (row k = column k of the original) so the advance can
	 * accumulate y += col_k * x[k] over contiguous, SIMD-friendly rows */
	double Vt[OPENPMX_STATE_MAX * OPENPMX_STATE_MAX];
	double Vinvt[OPENPMX_STATE_MAX * OPENPMX_STATE_MAX];

	/* per-mode constants, filled by eigen_decompose() */
	int degen[OPENPMX_STATE_MAX];				/* eigenvalue ~ 0: pure accumulation */
	double ieig[OPENPMX_STATE_MAX];				/* 1/(-eigval), 0 if degenerate */

	/* time-dependent cache, valid for one dt until re-decomposition.
	 * zeroed by advan_base_construct() so dtvalid starts at 0 */
	int dtvalid;
	double dtlast;
	double decay[OPENPMX_STATE_MAX];			/* exp(eigval*dt) */
	double step[OPENPMX_STATE_MAX];				/* (1-decay)/(-eigval), or dt if degenerate */
} ADVANCER_EIGEN;

static void advancer_eigen_info(const struct ADVANFUNCS* const advanfuncs,
                                FILE* f)
{
    fprintf(f, "advan model eigensystem (general linear)\n");
    fprintf(f, "advan nstate %i\n", advanfuncs->nstate);
}

static void advancer_eigen_construct(ADVAN* advan,
                                     const struct ADVANFUNCS* const advanfuncs)
{
	var self = container_of(advan, ADVANCER_EIGEN, advan);

	assert(advanfuncs->advan_size == sizeof(ADVANCER_EIGEN));
	advan_base_construct(&self->advan, advanfuncs);

	let n = advanfuncs->nstate;
	self->n = n;
	self->last_recalc_initcount = -1;
	
	self->sysmat_data = mallocvar(double, n * n);
	self->last_sysmat_data = mallocvar(double, n * n);

	/* Allocate the GSL eigen workspace once; reused at each decomposition */
	self->eigen_ws = gsl_eigen_nonsymmv_alloc(n);
	self->eval = gsl_vector_complex_alloc(n);
	self->evec = gsl_matrix_complex_alloc(n, n);    
	self->perm = gsl_permutation_alloc(n);
	self->eigvals = mallocvar(double, n);
	assert(self->eigen_ws);
	assert(self->eval);
	assert(self->evec);
	assert(self->perm);
	assert(self->eigvals);

	self->V_m = gsl_matrix_alloc(n, n);
	self->Vinv_m = gsl_matrix_alloc(n, n);
	assert(self->V_m);
	assert(self->Vinv_m);

	/* point to our sysmat so users can see it and update */
	advan->eigen_sysmat_data = self->sysmat_data;
}

static void advancer_eigen_destruct(ADVAN* advan)
{
	var self = container_of(advan, ADVANCER_EIGEN, advan);

	gsl_matrix_free(self->V_m);
	gsl_matrix_free(self->Vinv_m);

	free(self->eigvals);
	gsl_permutation_free(self->perm);
	gsl_vector_complex_free(self->eval);
	gsl_matrix_complex_free(self->evec);
	gsl_eigen_nonsymmv_free(self->eigen_ws);

	free(self->last_sysmat_data);
	free(self->sysmat_data);

	advan_base_destruct(advan);
}

/*------------------------------------------------------------------------
 * Eigendecompose the system matrix using GSL.
 *
 * Called once per parameter change. Fills eigvals, V, Vinv.
 * The system matrix S is factored as S = V * diag(eigvals) * Vinv.
 *
 * PK system matrices have real, negative eigenvalues (stable, distinct
 * modes). We assert this -- complex eigenvalues would indicate a
 * malformed PK model.
 *----------------------------------------------------------------------*/
__attribute__ ((hot))
static void eigen_decompose(ADVANCER_EIGEN* self, 
							const double sysmat_data[static OPENPMX_STATE_MAX * OPENPMX_STATE_MAX])
{
	/* Wrap in a GSL matrix view (row-major, which is GSL's native layout).
	 * gsl_eigen_nonsymmv destroys the input, so we work on a copy. */
	let n = self->n;
	double A_work[n * n];
	memcpy(A_work, sysmat_data, n * n * sizeof(double));
	var A_view = gsl_matrix_view_array(A_work, n, n);

	/* 2. Eigendecomposition via GSL */
	var eval = self->eval;
	var evec = self->evec;
	var status = gsl_eigen_nonsymmv(&A_view.matrix, eval, evec, self->eigen_ws);
	if (status != GSL_SUCCESS) {
		fprintf(stderr, "fatal: eigen ADVAN: gsl_eigen_nonsymmv failed "
				"(status=%d)\n", status);
		exit(EXIT_FAILURE);
	}

	/* Sort by ascending absolute value (smallest magnitude first).
	 * This gives a consistent ordering for the modes. */
	gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC); 

	/* 3. Extract real eigenvalues and eigenvectors */
	forcount(i, n) {
		let eval_i = gsl_vector_complex_get(eval, i);

		/* Verify eigenvalue is real (imaginary part negligible) */
		if (fabs(GSL_IMAG(eval_i)) > 1e-10) {
			fprintf(stderr, "fatal: eigen ADVAN: complex eigenvalue detected "
					"(lambda[%d] = %g + %gi). System matrix is not a "
					"valid PK model.\n", i, GSL_REAL(eval_i), GSL_IMAG(eval_i));
			exit(EXIT_FAILURE);
		}

		/* Verify eigenvalue is non-positive (stable system) */
		if (GSL_REAL(eval_i) > 1e-10) {
			fprintf(stderr, "fatal: eigen ADVAN: positive eigenvalue detected "
					"(lambda[%d] = %g). System is unstable.\n",
					i, GSL_REAL(eval_i));
			exit(EXIT_FAILURE);
		}

		self->eigvals[i] = GSL_REAL(eval_i);

		/* per-mode constants for the advance: one division here instead
		 * of one per mode per call */
		self->degen[i] = !(fabs(self->eigvals[i]) > 1e-15);
		self->ieig[i] = self->degen[i] ? 0. : -1./self->eigvals[i];

		/* Extract real parts of eigenvectors into V (row-major).
		 * GSL returns complex eigenvectors even for real eigenvalues;
		 * for PK systems the imaginary parts should be negligible. */
		forcount(j, n) {
			let v_ji = gsl_matrix_complex_get(evec, j, i);
			if (fabs(GSL_IMAG(v_ji)) > 1e-10) {
				fprintf(stderr, "fatal: eigen ADVAN: complex eigenvector "
						"component at (%d,%d). PK systems should have "
						"real eigenvectors.\n", j, i);
				exit(EXIT_FAILURE);
			}
			/* V[j, i] in row-major: row j, column i */
			gsl_matrix_set(self->V_m, j, i, GSL_REAL(v_ji));
		}
	}

	/* 4. Compute V^-1 via GSL LU decomposition */
	double Vcopy[n * n];
	var V_view = gsl_matrix_view_array(Vcopy, n, n);
	gsl_matrix_memcpy(&V_view.matrix, self->V_m);

	var perm = self->perm;
	int signum;
	status = gsl_linalg_LU_decomp(&V_view.matrix, perm, &signum);
	if (status != GSL_SUCCESS) {
		fprintf(stderr, "fatal: eigen ADVAN: LU decomposition of eigenvector "
				"matrix failed (status=%d).\n", status);
		exit(EXIT_FAILURE);
	}

	status = gsl_linalg_LU_invert(&V_view.matrix, perm, self->Vinv_m);
	if (status != GSL_SUCCESS) {
		fprintf(stderr, "fatal: eigen ADVAN: matrix inversion failed "
				"(status=%d).\n", status);
		exit(EXIT_FAILURE);
	}

	/* 5. transposed copies for the advance */
	forcount(i, n) {
		forcount(j, n) {
			self->Vt[j*n + i] = self->V_m->data[i*n + j];
			self->Vinvt[j*n + i] = self->Vinv_m->data[i*n + j];
		}
	}

#if 0
	/* 6. Verify reconstruction: S = V * diag(eigvals) * Vinv */
	forcount(i, n) {
		forcount(j, n) {
			double sum = 0.;
			for (int k = 0; k < n; k++) {
				/* V[i,k] * eigvals[k] * Vinv[k,j] (all row-major) */
				sum += self->V[i * n + k] * self->eigvals[k]
					 * self->Vinv[k * n + j];
			}
			if (fabs(sum - sysmat_data[i * n + j]) > 1e-6) {
				fprintf(stderr, "fatal: eigen ADVAN: reconstruction check "
						"failed at (%d,%d): got %g, expected %g\n",
						i, j, sum, sysmat_data[i * n + j]);
				exit(EXIT_FAILURE);
			}
		}
	}
#endif
}

/*------------------------------------------------------------------------
 * Advance interval: the core analytical state update.
 *
 * This implements the eigendecomposition equivalent of the STANPUMP
 * model() function, generalized to n compartments.
 *
 * Algorithm:
 *   1. Convert compartment state to modal coordinates:
 *        x_modal = Vinv * state
 *
 *   2. Project input (infusion rates) into modal coordinates:
 *        j_modal = Vinv * rates
 *
 *   3. Advance each mode analytically:
 *        x_modal_i(t+dt) = x_modal_i(t) * exp(eigval_i * dt)
 *                        + j_modal_i / (-eigval_i)
 *                          * (1 - exp(eigval_i * dt))
 *
 *      This is Jona's model() function:
 *        A(t+dt) = A(t) * decay + coef * rate * (1 - decay)
 *      applied independently to each mode.
 *
 *   4. Convert back to compartment coordinates:
 *        state = V * x_modal
 *
 * All matrices are row-major (GSL convention).
 *----------------------------------------------------------------------*/
__attribute__ ((hot))
static void advancer_eigen_advance_interval(ADVAN* advan,
											const IMODEL* const imodel,
											const RECORD* const record,
											double* const state,
											const POPPARAM* const popparam,
											const double endtime,
											const double* rates)
{
    var self = container_of(advan, ADVANCER_EIGEN, advan);

	(void)imodel;
	(void)record;
	(void)popparam;

	assert(advan->initcount > 0);

	/* Re-decompose only if system matrix changed. memcmp is
	 * SIMD-vectorized by libc and far faster than an element loop. */
	let n = self->n;
	if (self->last_recalc_initcount != advan->initcount) {
		if (memcmp(self->sysmat_data,
				   self->last_sysmat_data,
				   n * n * sizeof(double)) != 0) {
			eigen_decompose(self, self->sysmat_data);
			memcpy(self->last_sysmat_data,
				   self->sysmat_data,
				   n * n * sizeof(double));
			self->dtvalid = 0;
		}
		self->last_recalc_initcount = advan->initcount;
	}

	/* per-mode decay and step coefficients depend only on dt and the
	 * eigenvalues, so cache them for repeated dt */
	let dt = endtime - advan->time;
	assert(dt > 0.);
	if (!self->dtvalid || dt != self->dtlast) {
		for (int i = 0; i < n; i++) {
			if (__builtin_expect(self->degen[i], 0)) {
				/* limit of (1 - exp(eigval*dt)) / (-eigval) as eigval -> 0 */
				self->decay[i] = 1.;
				self->step[i] = dt;
			} else {
				let decay = exp(self->eigvals[i] * dt);
				self->decay[i] = decay;
				self->step[i] = (1. - decay) * self->ieig[i];
			}
		}
		self->dtlast = dt;
		self->dtvalid = 1;
	}

	/* modal advance on raw arrays: for n of 2-10 the GSL view/BLAS call
	 * overhead costs more than the arithmetic. Matrix-vector products are
	 * written column-wise (y[i] += Mt[k][i]*x[k]) so the inner loop over i
	 * vectorizes without -ffast-math and keeps the same summation order.
	 *   xm = decay .* (Vinv*state) + step .* (Vinv*rates)
	 *   state = V*xm */
	let Vinvt = self->Vinvt;
	let Vt = self->Vt;
	let decay = self->decay;
	let step = self->step;

	double xm[OPENPMX_STATE_MAX] = { 0. };
	double jm[OPENPMX_STATE_MAX] = { 0. };

	for (int k = 0; k < n; k++) {
		let col = Vinvt + k*n;
		let sk = state[k];
		for (int i = 0; i < n; i++)
			xm[i] += col[i]*sk;
	}

	/* skip the input projection when no infusion is running */
	int norates = 1;
	for (int k = 0; k < n; k++)
		if (rates[k] != 0.) {
			norates = 0;
			break;
		}
	if (norates) {
		for (int i = 0; i < n; i++)
			xm[i] *= decay[i];
	} else {
		for (int k = 0; k < n; k++) {
			let col = Vinvt + k*n;
			let rk = rates[k];
			for (int i = 0; i < n; i++)
				jm[i] += col[i]*rk;
		}
		for (int i = 0; i < n; i++)
			xm[i] = xm[i]*decay[i] + jm[i]*step[i];
	}

	/* back to compartment coordinates, in place (xm holds all inputs) */
	double out[OPENPMX_STATE_MAX] = { 0. };
	for (int k = 0; k < n; k++) {
		let col = Vt + k*n;
		let xk = xm[k];
		for (int i = 0; i < n; i++)
			out[i] += col[i]*xk;
	}
	for (int i = 0; i < n; i++)
		state[i] = out[i];
}

typedef struct {
	ADVANFUNCS advanfuncs;
} ADVANFUNCS_EIGEN;

ADVANFUNCS* pmx_advan_eigen(const DATACONFIG* const dataconfig,
							const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);
	assert(advanconfig->nstate > 0);
	
	var nstate = advanconfig->nstate;
	int err = advan_ensure(nstate >= 1, __func__, "nstate should be >= 1");
	if (err)
		return 0;

	let retinit = (ADVANFUNCS_EIGEN) {
		.advanfuncs = {
			.advan_size = sizeof(ADVANCER_EIGEN),
			.construct = advancer_eigen_construct,
			.destruct = advancer_eigen_destruct,
			.info = advancer_eigen_info,

			.reset = 0,
			.interval = advancer_eigen_advance_interval,

			.advanconfig = advanconfig,
			.recordinfo = recordinfo_init(dataconfig),
			.nstate = advanconfig->nstate,
		},
	};

	ADVANFUNCS* ret = malloc(sizeof(ADVANFUNCS_EIGEN));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANFUNCS_EIGEN));

	return ret;
}

/*------------------------------------------------------------------------
 *----------------------------------------------------------------------*/

static void advancer_eigen_threecomp_construct(ADVAN* advan,
											   const struct ADVANFUNCS* const advanfuncs)
{
	advancer_eigen_construct(advan,advanfuncs);
	
	/* hide the system matrix pointer that eigensolver exposes, because
	 * we do it ourselves. Any other use would be an error. */
	advan->eigen_sysmat_data = NULL;
}

static void advancer_eigen_threecomp_destruct(ADVAN* advan)
{
	advancer_eigen_destruct(advan);
}

static void advancer_eigen_threecomp_info(const struct ADVANFUNCS* const advanfuncs,
										  FILE* f)
{
    fprintf(f, "advan model eigensystem (three compartment)\n");
    fprintf(f, "advan nstate %i\n", advanfuncs->nstate);
}

typedef struct {
	ADVANFUNCS_EIGEN eigen;
	int offsetV1;
	int offsetV2;
	int offsetV3;
	int offsetCL;
	int offsetQ2;
	int offsetQ3;
} ADVANFUNCS_EIGEN_THREECOMP;

__attribute__ ((hot))
static void advancer_eigen_threecomp_advance_interval(ADVAN* advan,
													  const IMODEL* const imodel,
													  const RECORD* const record,
													  double* const state,
													  const POPPARAM* const popparam,
													  const double endtime,
													  const double* rates)
{
    var self = container_of(advan, ADVANCER_EIGEN, advan);
    
    /* do we need to update the eigensystem matrix? */
	if (self->last_recalc_initcount != advan->initcount) {
		var eigen_funcs   = container_of(advan->advanfuncs, ADVANFUNCS_EIGEN,           advanfuncs);
		var imodeloffsets = container_of(eigen_funcs,       ADVANFUNCS_EIGEN_THREECOMP, eigen); 
	 
		let V1 = *(const double*)(((char*)imodel) + imodeloffsets->offsetV1);
		let V2 = *(const double*)(((char*)imodel) + imodeloffsets->offsetV2);
		let V3 = *(const double*)(((char*)imodel) + imodeloffsets->offsetV3);
		let CL = *(const double*)(((char*)imodel) + imodeloffsets->offsetCL);
		let Q2 = *(const double*)(((char*)imodel) + imodeloffsets->offsetQ2);
		let Q3 = *(const double*)(((char*)imodel) + imodeloffsets->offsetQ3);
		
		/* define the eigensystem */
		let k10 = CL / V1;
		let k12 = Q2 / V1;
		let k21 = Q2 / V2;
		let k13 = Q3 / V1;
		let k31 = Q3 / V3;
		double sysmat_data[OPENPMX_STATE_MAX * OPENPMX_STATE_MAX] = { 
				-k10-k12-k13,	k21,	k31,
				k12,			-k21,	0.,
				k13, 			0., 	-k31 
		};
		
		/* update the system matrix for the eigen advance */
		let n = self->n;
		assert(n == advan->advanfuncs->nstate);
		memcpy(self->sysmat_data, sysmat_data, n * n * sizeof(double));

		/* resetting the last_recalc_initcount should be done in the 
		 * advancer_eigen_advance_interval() function, dont do it here */
	}

	/* pass the advance down to the eigensystem advancer */
	advancer_eigen_advance_interval(advan, imodel, record, state, popparam, endtime, rates);
}

/// This file also implements a three compartment mammilary model
/// via the underlying eigensystem solver. It is available as
/// `pmx_advan_eigen_threecomp()` of via openpmxtran as
/// `$ADVAN(eigen_threecomp)`.
ADVANFUNCS* pmx_advan_eigen_threecomp(const DATACONFIG* const dataconfig,
									  const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);
	
	var nstate = advanconfig->nstate;
	if (nstate == 0)
		nstate = 3;
	int err = advan_ensure(nstate >= 3, __func__, "nstate should be >= 3");
	if (err)
		return 0;

	let retinit = (ADVANFUNCS_EIGEN_THREECOMP) {
		.eigen = {
			.advanfuncs = {
				.advan_size = sizeof(ADVANCER_EIGEN),
				.construct = advancer_eigen_threecomp_construct,
				.destruct = advancer_eigen_threecomp_destruct,
				.info = advancer_eigen_threecomp_info,

				.reset = 0,
				.interval = advancer_eigen_threecomp_advance_interval,

				.advanconfig = advanconfig,
				.recordinfo = recordinfo_init(dataconfig),
				.nstate = nstate,
			},
		},
		.offsetV1 = structinfo_find_offset("V1", &advanconfig->imodelfields),
		.offsetV2 = structinfo_find_offset("V2", &advanconfig->imodelfields),
		.offsetV3 = structinfo_find_offset("V3", &advanconfig->imodelfields),
		.offsetCL = structinfo_find_offset("CL", &advanconfig->imodelfields),
		.offsetQ2 = structinfo_find_offset("Q2", &advanconfig->imodelfields),
		.offsetQ3 = structinfo_find_offset("Q3", &advanconfig->imodelfields),
	};
	err |= advan_ensure(retinit.offsetV1 >= 0, __func__, "could not find V1");
	err |= advan_ensure(retinit.offsetV2 >= 0, __func__, "could not find V2");
	err |= advan_ensure(retinit.offsetV3 >= 0, __func__, "could not find V3");
	err |= advan_ensure(retinit.offsetCL >= 0, __func__, "could not find CL");
	err |= advan_ensure(retinit.offsetQ2 >= 0, __func__, "could not find Q2");
	err |= advan_ensure(retinit.offsetQ3 >= 0, __func__, "could not find Q3");
	if (err)
		return 0;
	
	ADVANFUNCS* ret = malloc(sizeof(ADVANFUNCS_EIGEN_THREECOMP));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANFUNCS_EIGEN_THREECOMP));

	return ret;
}

/*------------------------------------------------------------------------
 *----------------------------------------------------------------------*/

static void advancer_eigen_twocomp_construct(ADVAN* advan,
											 const struct ADVANFUNCS* const advanfuncs)
{
	advancer_eigen_construct(advan,advanfuncs);
	
	/* hide the system matrix pointer that eigensolver exposes, because
	 * we do it ourselves. Any other use would be an error. */
	advan->eigen_sysmat_data = NULL;
}

static void advancer_eigen_twocomp_destruct(ADVAN* advan)
{
	advancer_eigen_destruct(advan);
}

static void advancer_eigen_twocomp_info(const struct ADVANFUNCS* const advanfuncs,
										FILE* f)
{
    fprintf(f, "advan model eigensystem (two compartment)\n");
    fprintf(f, "advan nstate %i\n", advanfuncs->nstate);
}

typedef struct {
	ADVANFUNCS_EIGEN eigen;
	int offsetV1;
	int offsetV2;
	int offsetCL; 
	int offsetQ2;
} ADVANFUNCS_EIGEN_TWOCOMP;

__attribute__ ((hot))
static void advancer_eigen_twocomp_advance_interval(ADVAN* advan,
													const IMODEL* const imodel,
													const RECORD* const record,
													double* const state,
													const POPPARAM* const popparam,
													const double endtime,
													const double* rates)
{
    var self = container_of(advan, ADVANCER_EIGEN, advan);
    
    /* do we need to update the eigensystem matrix? */
	if (self->last_recalc_initcount != advan->initcount) {
		var eigen_funcs   = container_of(advan->advanfuncs, ADVANFUNCS_EIGEN, advanfuncs);
		var imodeloffsets = container_of(eigen_funcs, ADVANFUNCS_EIGEN_TWOCOMP, eigen); 
	 
		let V1 = *(const double*)(((char*)imodel) + imodeloffsets->offsetV1);
		let V2 = *(const double*)(((char*)imodel) + imodeloffsets->offsetV2);
		let CL = *(const double*)(((char*)imodel) + imodeloffsets->offsetCL);
		let Q2 = *(const double*)(((char*)imodel) + imodeloffsets->offsetQ2);
		
		/* define the eigensystem */
		let k10 = CL / V1;
		let k12 = Q2 / V1;
		let k21 = Q2 / V2;
		double sysmat_data[OPENPMX_STATE_MAX * OPENPMX_STATE_MAX] = { 
			-k10-k12,  k21,
			k12,     -k21,
		};
		/* update the system matrix for the eigen advance */
		let n = self->n;
		assert(n == advan->advanfuncs->nstate);
		memcpy(self->sysmat_data, sysmat_data, n * n * sizeof(double));

		/* resetting the last_recalc_initcount should be done in the 
		 * advancer_eigen_advance_interval() function, dont do it here */
	}
	
	/* pass the advance down to the eigensystem advancer */
 	advancer_eigen_advance_interval(advan, imodel, record, state, popparam, endtime, rates);
}

/// This file also implements a two compartment mammilary model via the
/// underlying eigensystem solver. It is available as
/// `pmx_advan_eigen_twocomp()` of via openpmxtran as
/// `$ADVAN(eigen_twocomp)`.
ADVANFUNCS* pmx_advan_eigen_twocomp(const DATACONFIG* const dataconfig,
									const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);
	
	var nstate = advanconfig->nstate;
	if (nstate == 0)
		nstate = 2;
	int err = advan_ensure(nstate >= 2, __func__, "nstate should be >= 2");
	if (err)
		return 0;
	
	let retinit = (ADVANFUNCS_EIGEN_TWOCOMP) {
		.eigen = {
			.advanfuncs = {
				.advan_size = sizeof(ADVANCER_EIGEN),
				.construct = advancer_eigen_twocomp_construct,
				.destruct = advancer_eigen_twocomp_destruct,
				.info = advancer_eigen_twocomp_info,

				.reset = 0,
				.interval = advancer_eigen_twocomp_advance_interval,

				.advanconfig = advanconfig,
				.recordinfo = recordinfo_init(dataconfig),
				.nstate = nstate,
			},
		},
		.offsetV1 = structinfo_find_offset("V1", &advanconfig->imodelfields),
		.offsetV2 = structinfo_find_offset("V2", &advanconfig->imodelfields),
		.offsetCL = structinfo_find_offset("CL", &advanconfig->imodelfields),
		.offsetQ2 = structinfo_find_offset("Q2", &advanconfig->imodelfields),
	};
	err |= advan_ensure(retinit.offsetV1 >= 0, __func__, "could not find V1");
	err |= advan_ensure(retinit.offsetV2 >= 0, __func__, "could not find V2");
	err |= advan_ensure(retinit.offsetCL >= 0, __func__, "could not find CL");
	err |= advan_ensure(retinit.offsetQ2 >= 0, __func__, "could not find Q2");
	if (err)
		return 0;
	
	ADVANFUNCS* ret = malloc(sizeof(ADVANFUNCS_EIGEN_TWOCOMP));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANFUNCS_EIGEN_TWOCOMP));

	return ret;
}

/*------------------------------------------------------------------------
 *----------------------------------------------------------------------*/

static void advancer_eigen_onecomp_depot_construct(ADVAN* advan,
													const struct ADVANFUNCS* const advanfuncs)
{
	advancer_eigen_construct(advan,advanfuncs);
	
	/* hide the system matrix pointer that eigensolver exposes, because
	 * we do it ourselves. Any other use would be an error. */
	advan->eigen_sysmat_data = NULL;
}

static void advancer_eigen_onecomp_depot_destruct(ADVAN* advan)
{
	advancer_eigen_destruct(advan);
}

static void advancer_eigen_onecomp_depot_info(const struct ADVANFUNCS* const advanfuncs,
											   FILE* f)
{
    fprintf(f, "advan model eigensystem (one compartment with depot)\n");
    fprintf(f, "advan nstate %i\n", advanfuncs->nstate);
}

typedef struct {
	ADVANFUNCS_EIGEN eigen;
	int offsetV;
	int offsetCL; 
	int offsetKA;
} ADVANFUNCS_EIGEN_ONECOMP_DEPOT;

__attribute__ ((hot))
static void advancer_eigen_onecomp_depot_advance_interval(ADVAN* advan,
														   const IMODEL* const imodel,
														   const RECORD* const record,
														   double* const state,
														   const POPPARAM* const popparam,
														   const double endtime,
														   const double* rates)
{
    var self = container_of(advan, ADVANCER_EIGEN, advan);
    
    /* do we need to update the eigensystem matrix? */
	if (self->last_recalc_initcount != advan->initcount) {
		var eigen_funcs   = container_of(advan->advanfuncs, ADVANFUNCS_EIGEN, advanfuncs);
		var imodeloffsets = container_of(eigen_funcs, ADVANFUNCS_EIGEN_ONECOMP_DEPOT, eigen); 
	 
		let V = *(const double*)(((char*)imodel) + imodeloffsets->offsetV);
		let CL = *(const double*)(((char*)imodel) + imodeloffsets->offsetCL);
		let KA = *(const double*)(((char*)imodel) + imodeloffsets->offsetKA);
		
		/* define the eigensystem */
		let KE = CL / V;
		double sysmat_data[OPENPMX_STATE_MAX * OPENPMX_STATE_MAX] = { 
			-KA,	0.,
			KA,		-KE,
		};
		/* update the system matrix for the eigen advance */
		let n = self->n;
		assert(n == advan->advanfuncs->nstate);
		memcpy(self->sysmat_data, sysmat_data, n * n * sizeof(double));

		/* resetting the last_recalc_initcount should be done in the 
		 * advancer_eigen_advance_interval() function, dont do it here */
	}
	
	/* pass the advance down to the eigensystem advancer */
 	advancer_eigen_advance_interval(advan, imodel, record, state, popparam, endtime, rates);
}

/// This file also implements a one compartment model with depot
/// via the underlying eigensystem solver. It is available as
/// `pmx_advan_eigen_onecomp_depot()` or via openpmxtran as
/// `$ADVAN(eigen_onecomp_depot)`.
ADVANFUNCS* pmx_advan_eigen_onecomp_depot(const DATACONFIG* const dataconfig,
										  const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);
	
	var nstate = advanconfig->nstate;
	if (nstate == 0)
		nstate = 2;
	int err = advan_ensure(nstate >= 2, __func__, "nstate should be >= 2");
	if (err)
		return 0;

	let retinit = (ADVANFUNCS_EIGEN_ONECOMP_DEPOT) {
		.eigen = {
			.advanfuncs = {
				.advan_size = sizeof(ADVANCER_EIGEN),
				.construct = advancer_eigen_onecomp_depot_construct,
				.destruct = advancer_eigen_onecomp_depot_destruct,
				.info = advancer_eigen_onecomp_depot_info,

				.reset = 0,
				.interval = advancer_eigen_onecomp_depot_advance_interval,

				.advanconfig = advanconfig,
				.recordinfo = recordinfo_init(dataconfig),
				.nstate = nstate,
			},
		},
		.offsetV = structinfo_find_offset("V", &advanconfig->imodelfields),
		.offsetCL = structinfo_find_offset("CL", &advanconfig->imodelfields),
		.offsetKA = structinfo_find_offset("KA", &advanconfig->imodelfields),
	};
	err |= advan_ensure(retinit.offsetV >= 0, __func__, "could not find V");
	err |= advan_ensure(retinit.offsetCL >= 0, __func__, "could not find CL");
	err |= advan_ensure(retinit.offsetKA >= 0, __func__, "could not find KA");
	if (err)
		return 0;
	
	ADVANFUNCS* ret = malloc(sizeof(ADVANFUNCS_EIGEN_ONECOMP_DEPOT));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANFUNCS_EIGEN_ONECOMP_DEPOT));

	return ret;
}


/*------------------------------------------------------------------------
 *----------------------------------------------------------------------*/

static void advancer_eigen_twocomp_depot_construct(ADVAN* advan,
													const struct ADVANFUNCS* const advanfuncs)
{
	advancer_eigen_construct(advan,advanfuncs);
	
	/* hide the system matrix pointer that eigensolver exposes, because
	 * we do it ourselves. Any other use would be an error. */
	advan->eigen_sysmat_data = NULL;
}

static void advancer_eigen_twocomp_depot_destruct(ADVAN* advan)
{
	advancer_eigen_destruct(advan);
}

static void advancer_eigen_twocomp_depot_info(const struct ADVANFUNCS* const advanfuncs,
											   FILE* f)
{
    fprintf(f, "advan model eigensystem (two compartment with depot)\n");
    fprintf(f, "advan nstate %i\n", advanfuncs->nstate);
}

typedef struct {
	ADVANFUNCS_EIGEN eigen;
	int offsetV1;
	int offsetV2;
	int offsetCL; 
	int offsetQ2; 
	int offsetKA;
} ADVANFUNCS_EIGEN_TWOCOMP_DEPOT;

__attribute__ ((hot))
static void advancer_eigen_twocomp_depot_advance_interval(ADVAN* advan,
														   const IMODEL* const imodel,
														   const RECORD* const record,
														   double* const state,
														   const POPPARAM* const popparam,
														   const double endtime,
														   const double* rates)
{
    var self = container_of(advan, ADVANCER_EIGEN, advan);
    
    /* do we need to update the eigensystem matrix? */
	if (self->last_recalc_initcount != advan->initcount) {
		var eigen_funcs   = container_of(advan->advanfuncs, ADVANFUNCS_EIGEN, advanfuncs);
		var imodeloffsets = container_of(eigen_funcs, ADVANFUNCS_EIGEN_TWOCOMP_DEPOT, eigen); 
	 
		let V1 = *(const double*)(((char*)imodel) + imodeloffsets->offsetV1);
		let V2 = *(const double*)(((char*)imodel) + imodeloffsets->offsetV2);
		let CL = *(const double*)(((char*)imodel) + imodeloffsets->offsetCL);
		let Q2 = *(const double*)(((char*)imodel) + imodeloffsets->offsetQ2);
		let KA = *(const double*)(((char*)imodel) + imodeloffsets->offsetKA);
		
		/* define the eigensystem */
		double K10=CL/V1;
		double K12=Q2/V1;
		double K21=Q2/V2;
		double sysmat_data[OPENPMX_STATE_MAX * OPENPMX_STATE_MAX] = { 
			-KA,		0., 		0.,
			KA,		-(K10+K12),	K21,
			0.,		K12, 		-K21 };

		/* update the system matrix for the eigen advance */
		let n = self->n;
		assert(n == advan->advanfuncs->nstate);
		memcpy(self->sysmat_data, sysmat_data, n * n * sizeof(double));

		/* resetting the last_recalc_initcount should be done in the 
		 * advancer_eigen_advance_interval() function, dont do it here */
	}
	
	/* pass the advance down to the eigensystem advancer */
 	advancer_eigen_advance_interval(advan, imodel, record, state, popparam, endtime, rates);
}

/// This file also implements a two compartment model with depot
/// via the underlying eigensystem solver. It is available as
/// `pmx_advan_eigen_twocomp_depot()` or via openpmxtran as
/// `$ADVAN(eigen_twocomp_depot)`.
ADVANFUNCS* pmx_advan_eigen_twocomp_depot(const DATACONFIG* const dataconfig,
										  const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);

	var nstate = advanconfig->nstate;
	if (nstate == 0)
		nstate = 3;
	int err = advan_ensure(nstate >= 3, __func__, "nstate should be >= 3");
	if (err)
		return 0;

	let retinit = (ADVANFUNCS_EIGEN_TWOCOMP_DEPOT) {
		.eigen = {
			.advanfuncs = {
				.advan_size = sizeof(ADVANCER_EIGEN),
				.construct = advancer_eigen_twocomp_depot_construct,
				.destruct = advancer_eigen_twocomp_depot_destruct,
				.info = advancer_eigen_twocomp_depot_info,

				.reset = 0,
				.interval = advancer_eigen_twocomp_depot_advance_interval,

				.advanconfig = advanconfig,
				.recordinfo = recordinfo_init(dataconfig),
				.nstate = nstate,
			},
		},
		.offsetV1 = structinfo_find_offset("V1", &advanconfig->imodelfields),
		.offsetV2 = structinfo_find_offset("V2", &advanconfig->imodelfields),
		.offsetCL = structinfo_find_offset("CL", &advanconfig->imodelfields),
		.offsetQ2 = structinfo_find_offset("Q2", &advanconfig->imodelfields),
		.offsetKA = structinfo_find_offset("KA", &advanconfig->imodelfields),
	};
	err |= advan_ensure(retinit.offsetV1 >= 0, __func__, "could not find V1");
	err |= advan_ensure(retinit.offsetV2 >= 0, __func__, "could not find V2");
	err |= advan_ensure(retinit.offsetCL >= 0, __func__, "could not find CL");
	err |= advan_ensure(retinit.offsetQ2 >= 0, __func__, "could not find Q2");
	err |= advan_ensure(retinit.offsetKA >= 0, __func__, "could not find KA");
	if (err)
		return 0;
	
	ADVANFUNCS* ret = malloc(sizeof(ADVANFUNCS_EIGEN_TWOCOMP_DEPOT));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANFUNCS_EIGEN_TWOCOMP_DEPOT));

	return ret;
}



