/*
 * Eigendecomposition ADVAN for OpenPMX
 *
 * This advancer uses eigendecomposition to obtain exact analytical
 * solutions for any linear compartmental model, regardless of the
 * number of compartments or topology.
 *
 * Instead of hardcoding the analytical solution for a specific model
 * structure (like onecomp.c or threecomp.c), or using numerical ODE
 * integration (like diffeqn_libgsl.c), this advancer:
 *
 *   1. User must supply system matrix
 *   2. Eigendecomposes the matrix: A = V * diag(eigvals) * V^-1
 *   3. Advances state analytically in modal coordinates
 *
 * This is equivalent to NONMEM's ADVAN5 (linear general compartmental
 * model). The technique is described in:
 *   - Joachim J, "TCI Eigendecomposition" (2026)
 *   - Minto C, Schnider T, "From Eigenvalues to Plasma Coefficients" (2026)
 *
 * The modal state update for constant input over interval dt:
 *
 *   x_modal(t+dt) = x_modal(t) * exp(eigval*dt)
 *                  + (V^-1 * J) / (-eigval) * (1 - exp(eigval*dt))
 *
 * where eigval < 0 for stable systems, V is the eigenvector matrix,
 * and J is the input vector (infusion rates into each compartment).
 *
 * Uses GSL for eigendecomposition (gsl_eigen_nonsymmv) and matrix
 * inversion (LU decomposition).
 * Link with: -lgsl -lgslcblas
 */

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
#include <gsl/gsl_blas.h>

/*------------------------------------------------------------------------
 * sysmat is stored in row-major order (C convention, GSL convention):
 *   sysmat[i * n + j] = element at row i, column j
 *
 * The matrix operates on amounts (consistent with OpenPMX state).
 * Infusion rates are handled separately by the advancer.
 *----------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 * Advancer state: extends ADVAN with eigen-specific cached data.
 *----------------------------------------------------------------------*/
typedef struct {
	ADVAN advan;

	/* system matrix */
	/* user will fill this in in $IMODEL() call, via SYSMAT() macro */
	double sysmat_data[OPENPMX_STATE_MAX * OPENPMX_STATE_MAX];
	/* cached sysmat to avoid unnecessary recalculation */
	double last_sysmat_data[OPENPMX_STATE_MAX * OPENPMX_STATE_MAX];
	int last_recalc_initcount;
	
	/* Cached eigendecomposition results */
	int n;                             					/* matrix dimension */
	double eigvals[OPENPMX_STATE_MAX];         			/* eigenvalues (real parts) */
	double V[OPENPMX_STATE_MAX * OPENPMX_STATE_MAX];  	/* right eigenvector matrix */
	double Vinv[OPENPMX_STATE_MAX * OPENPMX_STATE_MAX];	/* inverse of V */

	/* GSL workspaces allocated once in constructor, reused */
	gsl_eigen_nonsymmv_workspace* eigen_ws;
	gsl_vector_complex* eval;
	gsl_matrix_complex* evec;
	gsl_permutation* perm;    
} ADVANCER_EIGEN;

/*------------------------------------------------------------------------
 * Info
 *----------------------------------------------------------------------*/
static void advancer_eigen_info(const struct ADVANFUNCS* const advanfuncs,
                                FILE* f)
{
    fprintf(f, "advan model eigensystem (general linear)\n");
    fprintf(f, "advan nstate %i\n", advanfuncs->nstate);
}

/*------------------------------------------------------------------------
 * Constructor: allocate GSL workspaces
 *----------------------------------------------------------------------*/
static void advancer_eigen_construct(ADVAN* advan,
                                     const struct ADVANFUNCS* const advanfuncs)
{
	var self = container_of(advan, ADVANCER_EIGEN, advan);

	assert(advanfuncs->advan_size == sizeof(ADVANCER_EIGEN));
	advan_base_construct(&self->advan, advanfuncs);

	forcount(i, OPENPMX_STATE_MAX * OPENPMX_STATE_MAX) {
		self->sysmat_data[i] = NAN;
		self->last_sysmat_data[i] = NAN;
	}
	self->n = advanfuncs->nstate;
	self->last_recalc_initcount = -1;

	/* Allocate the GSL eigen workspace once; reused at each decomposition */
	self->eigen_ws = gsl_eigen_nonsymmv_alloc(self->n);
	self->eval = gsl_vector_complex_alloc(self->n);
	self->evec = gsl_matrix_complex_alloc(self->n, self->n);    
	self->perm = gsl_permutation_alloc(self->n);
	assert(self->eigen_ws);
	assert(self->eval);
	assert(self->evec);
	assert(self->perm);

	/* point to our sysmat so users can see it and update */
	advan->eigen_sysmat_data = self->sysmat_data;
}

/*------------------------------------------------------------------------
 * Destructor: free GSL workspaces
 *----------------------------------------------------------------------*/
static void advancer_eigen_destruct(ADVAN* advan)
{
	var self = container_of(advan, ADVANCER_EIGEN, advan);

	gsl_permutation_free(self->perm);
	gsl_vector_complex_free(self->eval);
	gsl_matrix_complex_free(self->evec);
	gsl_eigen_nonsymmv_free(self->eigen_ws);

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
	double A_work[OPENPMX_STATE_MAX * OPENPMX_STATE_MAX];
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
			self->V[j * n + i] = GSL_REAL(v_ji);
		}
	}

	/* 4. Compute V^-1 via GSL LU decomposition */
	double Vcopy[OPENPMX_STATE_MAX * OPENPMX_STATE_MAX];
	memcpy(Vcopy, self->V, n * n * sizeof(double));
	var V_view = gsl_matrix_view_array(Vcopy, n, n);

	var perm = self->perm;
	int signum;
	status = gsl_linalg_LU_decomp(&V_view.matrix, perm, &signum);
	if (status != GSL_SUCCESS) {
		fprintf(stderr, "fatal: eigen ADVAN: LU decomposition of eigenvector "
				"matrix failed (status=%d).\n", status);
		exit(EXIT_FAILURE);
	}

	var Vinv_view = gsl_matrix_view_array(self->Vinv, n, n);
	status = gsl_linalg_LU_invert(&V_view.matrix, perm, &Vinv_view.matrix);
	if (status != GSL_SUCCESS) {
		fprintf(stderr, "fatal: eigen ADVAN: matrix inversion failed "
				"(status=%d).\n", status);
		exit(EXIT_FAILURE);
	}

#if 0
	/* 5. Verify reconstruction: S = V * diag(eigvals) * Vinv */
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
	 * SIMD-vectorized by libc and far faster than the previous scalar
	 * element loop. */
	let n = self->n;
	if (self->last_recalc_initcount != advan->initcount) {
		if (memcmp(self->sysmat_data,
				   self->last_sysmat_data,
				   n * n * sizeof(double)) != 0) {
			eigen_decompose(self, self->sysmat_data);
			memcpy(self->last_sysmat_data,
				   self->sysmat_data,
				   n * n * sizeof(double));
		}
		self->last_recalc_initcount = advan->initcount;
	}

	/* Wrap cached matrices and working vectors as GSL views.
	 * gsl_matrix_view_array / gsl_vector_view_array are zero-cost: they
	 * only fill a small header struct pointing at the existing arrays. */
	var Vinv_m = gsl_matrix_view_array(self->Vinv, n, n);
	var V_m    = gsl_matrix_view_array(self->V, n, n);

	/* Steps 1+2: project state and rates into modal coordinates via Vinv.
	 *
	 * x_modal = Vinv * state
	 * j_modal = Vinv * rates
	 *
	 * cblas_dgemv with CblasRowMajor / CblasNoTrans computes:
	 *   y = alpha * A * x + beta * y
	 * where A is (n x n) row-major, x and y are length-n vectors.
	 * alpha=1, beta=0 gives a plain matrix-vector product.
	 *
	 * The two calls share the same matrix Vinv; BLAS will use optimised
	 * SIMD kernels (SSE2/AVX on x86, NEON on ARM) automatically.
	 */
	double x_modal[OPENPMX_STATE_MAX];
	double j_modal[OPENPMX_STATE_MAX];
	var x_modal_v = gsl_vector_view_array(x_modal, n);
	var j_modal_v = gsl_vector_view_array(j_modal, n);
	var state_v   = gsl_vector_view_array(state, n);
	/* rates is const double*; gsl_vector_const_view_array avoids the cast */
	var rates_v   = gsl_vector_const_view_array(rates, n);

	gsl_blas_dgemv(CblasNoTrans, 1.0, &Vinv_m.matrix,
	               &state_v.vector, 0.0, &x_modal_v.vector);   /* x_modal = Vinv * state */
	gsl_blas_dgemv(CblasNoTrans, 1.0, &Vinv_m.matrix,
	               &rates_v.vector, 0.0, &j_modal_v.vector);   /* j_modal = Vinv * rates */

	/* Step 3: advance each mode analytically.
	 *
	 * For each eigenvalue lambda_i (negative for stable systems):
	 *   decay = exp(lambda_i * dt)
	 *   x_new = x_old * decay + j_modal * (1 - decay) / (-lambda_i)
	 *
	 * The (1 - decay)/(-lambda) term is the step response coefficient.
	 * This is exactly STANPUMP's p_coef mechanism generalized to n
	 * dimensions -- see Minto's equation (57).
	 */
	let dt = endtime - advan->time;
	assert(dt > 0.);
	double x_modal_new[OPENPMX_STATE_MAX];
	for (int i = 0; i < n; i++) {
		let eigval = self->eigvals[i];
		let decay  = exp(eigval * dt);

		if (fabs(eigval) > 1e-15) {
			/* Normal case: mode decays and accumulates */
			let step_coef = (1.0 - decay) / (-eigval);
			x_modal_new[i] = x_modal[i] * decay + j_modal[i] * step_coef;
		} else {
			/* Degenerate case: eigenvalue ~ 0, pure accumulation.
			 * Limit as eigval -> 0: (1 - exp(eigval*dt)) / (-eigval) -> dt */
			x_modal_new[i] = x_modal[i] + j_modal[i] * dt;
		}
	}

	/* Step 4: transform back to compartment coordinates.
	 *
	 * state = V * x_modal_new
	 *
	 * state_v still points at the `state` array so the result is written
	 * back in-place — no extra copy needed.
	 */
	var x_modal_new_v = gsl_vector_view_array(x_modal_new, n);
	gsl_blas_dgemv(CblasNoTrans, 1.0, &V_m.matrix,
	               &x_modal_new_v.vector, 0.0, &state_v.vector); /* state = V * x_modal_new */
}

/*------------------------------------------------------------------------
 * Public constructor: pmx_advan_eigen
 *
 * This advancer replaces both the hardcoded analytical ADVANs and the
 * ODE solver for any linear model. It is exact (no step size error)
 * and general (any number of compartments, any topology).
 *----------------------------------------------------------------------*/

typedef struct {
	ADVANFUNCS advanfuncs;
} ADVANFUNCS_EIGEN;

ADVANFUNCS* pmx_advan_eigen(const DATACONFIG* const dataconfig,
							const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);
	assert(advanconfig->nstate > 0);
	assert(advanconfig->nstate <= OPENPMX_STATE_MAX);

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
 * Public constructor: pmx_advan_eigen_threecomp
 *----------------------------------------------------------------------*/

static void advancer_eigen_threecomp_construct(ADVAN* advan,
											   const struct ADVANFUNCS* const advanfuncs)
{
	advancer_eigen_construct(advan,advanfuncs);
	
	/* hide the system matrix pointer that exposes, because we do it 
	 * ourselves. Any other use would be an error. */
	advan->eigen_sysmat_data = 0;
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
 	advancer_eigen_advance_interval(advan, imodel, record, state, popparam, endtime, rates);
}

ADVANFUNCS* pmx_advan_eigen_threecomp(const DATACONFIG* const dataconfig,
									  const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);
	assert(advanconfig->nstate == 0 || advanconfig->nstate == 3);

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
				.nstate = 3, /* dont allow user to set */
			},
		},
		.offsetV1 = structinfo_find_offset("V1", &advanconfig->imodelfields),
		.offsetV2 = structinfo_find_offset("V2", &advanconfig->imodelfields),
		.offsetV3 = structinfo_find_offset("V3", &advanconfig->imodelfields),
		.offsetCL = structinfo_find_offset("CL", &advanconfig->imodelfields),
		.offsetQ2 = structinfo_find_offset("Q2", &advanconfig->imodelfields),
		.offsetQ3 = structinfo_find_offset("Q3", &advanconfig->imodelfields),
	};
	advan_ensure(retinit.offsetV1 >= 0, __func__, "could not find V1");
	advan_ensure(retinit.offsetV2 >= 0, __func__, "could not find V2");
	advan_ensure(retinit.offsetV3 >= 0, __func__, "could not find V3");
	advan_ensure(retinit.offsetCL >= 0, __func__, "could not find CL");
	advan_ensure(retinit.offsetQ2 >= 0, __func__, "could not find Q2");
	advan_ensure(retinit.offsetQ3 >= 0, __func__, "could not find Q3");
	
	ADVANFUNCS* ret = malloc(sizeof(ADVANFUNCS_EIGEN_THREECOMP));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANFUNCS_EIGEN_THREECOMP));

	return ret;
}

/*------------------------------------------------------------------------
 * Public constructor: pmx_advan_eigen_twocomp
 *----------------------------------------------------------------------*/

static void advancer_eigen_twocomp_construct(ADVAN* advan,
											   const struct ADVANFUNCS* const advanfuncs)
{
	advancer_eigen_construct(advan,advanfuncs);
	
	/* hide the system matrix pointer that exposes, because we do it 
	 * ourselves. Any other use would be an error. */
	advan->eigen_sysmat_data = 0;
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
 	advancer_eigen_advance_interval(advan, imodel, record, state, popparam, endtime, rates);
}

ADVANFUNCS* pmx_advan_eigen_twocomp(const DATACONFIG* const dataconfig,
									const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);
	assert(advanconfig->nstate == 0 || advanconfig->nstate == 2);

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
				.nstate = 2, /* dont allow user to set */
			},
		},
		.offsetV1 = structinfo_find_offset("V1", &advanconfig->imodelfields),
		.offsetV2 = structinfo_find_offset("V2", &advanconfig->imodelfields),
		.offsetCL = structinfo_find_offset("CL", &advanconfig->imodelfields),
		.offsetQ2 = structinfo_find_offset("Q2", &advanconfig->imodelfields),
	};
	advan_ensure(retinit.offsetV1 >= 0, __func__, "could not find V1");
	advan_ensure(retinit.offsetV2 >= 0, __func__, "could not find V2");
	advan_ensure(retinit.offsetCL >= 0, __func__, "could not find CL");
	advan_ensure(retinit.offsetQ2 >= 0, __func__, "could not find Q2");
	
	ADVANFUNCS* ret = malloc(sizeof(ADVANFUNCS_EIGEN_TWOCOMP));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANFUNCS_EIGEN_TWOCOMP));

	return ret;
}


