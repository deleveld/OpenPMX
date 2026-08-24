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
#include "utils/various.h"

#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_eigen.h>

/* TODO: what should be done about number of evaluations, does
 * covariance count into PMX object?
 * 
 * How can I return information to the user?
 * 
 * Should I adjust H per dimension to have a specific delta objfn?

Need to make a diagnostics script.
	For diagnostics see:
	https://chatgpt.com/share/6a8ab577-0768-83ed-bfee-af831ae642c9

 * */

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
	popmodel->result.nsig = 0.;
}

static void decode_popmodel_limit(ENCODE* const encode, const gsl_matrix* const Hinv, const int j, const double z) 
{
	/* Gemini  says this */
//	let variance = 2. * gsl_matrix_get(Hinv, j, j);
	/* Claude and Chat GPT says this */
	let variance = 4. * gsl_matrix_get(Hinv, j, j);	
	let SD = sqrt(variance);

	typeof(encode->offset) offset = { 0 };

	offset[j] = z * SD;
	encode_untransform(encode, offset);
}

static FILE* estimate_results_fopen(const char* name,
									const char* ext)
{
	char fname[PATH_MAX];
	snprintf(fname, sizeof(fname), "%s%s", name, ext);
	return fopen(fname, "w");
}

void pmx_covariance(OPENPMX* const source, COVARIANCECONFIG* const userargs)
{
	gsl_matrix *H = 0;
	gsl_matrix *Hinv = 0;
	double* warm_etas = 0;
	
	/* Work on a copy of the source object */
	var pmx = pmx_copy(source);
	pmx_ensure_state(&pmx);
	var pstate = pmx.state;
	
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

	/* open the results file */
	FILE* outstream = 0;
	if (source->filename)
		outstream = estimate_results_fopen(source->filename, OPENPMX_COVFILE);
	info(outstream, "OpenPMX %i.%i.%i hash %s\n", OPENPMX_VERSION_MAJOR, OPENPMX_VERSION_MINOR, OPENPMX_VERSION_RELEASE, OPENPMX_GITHASH);

	/* keep a copy of the eta values so we can have a warm start */
	let idata = &pstate->idata;
	let nindivid = idata->nindivid;
	warm_etas = idata_etas_copy_alloc(idata);

	ERRCTX errctx = { 0 };
	let popmodel = popmodel_init(pmx.theta, pmx.omega, pmx.sigma, &errctx);
	if (errctx.len) 
		fatal(0, "cov: %s: %s", __func__, errctx.errmsg);

	var encode = encode_init(&popmodel);
	encode_transform(&encode, &popmodel);
	let nparam = encode.nparam;

	info(outstream, "cov nparam %i step_size %g\n", nparam, args.step_size);

	if (args.alpha >= 1.) {
		info(outstream, "cov: error: alpha >= 1, covariance not calculated\n");
		goto failed;
	}
	if (nparam == 0) {
		info(outstream, "cov: error: no estimated parameters, covariance not calculated\n");
		goto failed;
	}
	if (nindivid < nparam) {
		info(outstream, "cov: error: nindivid < nparam, covariance not calculated\n");
		goto failed;
	}

	var J = gsl_matrix_alloc(nindivid, nparam);
	var iobjfn_plus = mallocvar(double, nindivid);
	var iobjfn_minus = mallocvar(double, nindivid);

	/* Single Pass: Compute OPG Jacobian along canonical coordinate axes
	 * also known as the BHHH or OPG/Fisher information estimate */
	forcount(j, nparam) {
		let h = args.step_size;
		typeof(encode.offset) offset = { 0 };

		/* per-individual eta reset */
		idata_etas_set(idata, warm_etas);

		/* Evaluate f(+h_j) in transformed space */
		offset[j] = h;
		encode_untransform(&encode, offset);
		encode_evaluate(&encode, idata, pstate->advanfuncs, &options);
		forcount(k, nindivid)
			iobjfn_plus[k] = idata->individ[k].iobjfn;

		/* show progress */
		struct timespec now;
		clock_gettime(CLOCK_MONOTONIC, &now);
		var t = timespec_time_difference(&begin, &now) / 1000.;
		popmodel_eval_information(&encode.popmodel, t, 0, false, outstream, 0);

		/* per-individual eta reset */
		idata_etas_set(idata, warm_etas);

		/* Evaluate f(-h_j) in transformed space */
		offset[j] = -h;
		encode_untransform(&encode, offset);
		encode_evaluate(&encode, idata, pstate->advanfuncs, &options);
		forcount(k, nindivid) 
			iobjfn_minus[k] = idata->individ[k].iobjfn;;

		/* show progress */
		clock_gettime(CLOCK_MONOTONIC, &now);
		t = timespec_time_difference(&begin, &now) / 1000.;
		popmodel_eval_information(&encode.popmodel, t, 0, false, outstream, 0);

		/* Compute central difference gradient for column j */
		forcount(k, nindivid) {
			let v = (iobjfn_plus[k] - iobjfn_minus[k]) / (2.0 * h);
			gsl_matrix_set(J, k, j, v);
		}
	}
	free(iobjfn_plus);
	free(iobjfn_minus);

	/* Form approximate Hessian H = J^T * J */
	/* also known as BHHH method or outer-product-gradient approximation */
	H = gsl_matrix_alloc(nparam, nparam);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, J, J, 0.0, H);
	gsl_matrix_free(J);

	/* invert approximate Hessian */
	Hinv = gsl_matrix_alloc(nparam, nparam);
	gsl_matrix_memcpy(Hinv, H);
	if (gsl_linalg_cholesky_decomp(Hinv) != GSL_SUCCESS ||
		gsl_linalg_cholesky_invert(Hinv) != GSL_SUCCESS) {
		info(outstream, "cov: error: singular or ill-conditioned OPG matrix\n");
		goto failed;
	}

	// Critical p-value threshold
	// Split alpha in half to account for both the left and right tails
	let z = gsl_cdf_ugaussian_Qinv(args.alpha / 2.0);
	info(outstream, "cov alpha %g Z %g\n", args.alpha, z);

	info(outstream, OPENPMX_SFORMAT " " OPENPMX_SFORMAT OPENPMX_SFORMAT "  " OPENPMX_SFORMAT "\n",
		"PARAMETER", "LOWER", "VALUE", "UPPER ");

	/* This is dependant on the order of the parameter offsets in the
	 * ENCODE object. It would perhaps be better if there was a way
	 * to indicate the relationship between the transformed and
	 * untransformed spaces */
	char name[64], lower[64], upper[64];
	var j = 0;
	forcount(i, popmodel.ntheta) {
		if (popmodel.thetaestim[i] != FIXED) {
			decode_popmodel_limit(&encode, Hinv, j, z);
			snprintf(upper, sizeof(upper), OPENPMX_FFORMAT, encode.popmodel.theta[i]);

			decode_popmodel_limit(&encode, Hinv, j, -z); 
			snprintf(lower, sizeof(lower), OPENPMX_FFORMAT, encode.popmodel.theta[i]);

			let value = popmodel.theta[i];
			let ind = (source->data._offset1) ? (i + 1) : (i);
			snprintf(name, sizeof(name), "THETA(%i)", ind); 
			info(outstream, OPENPMX_SFORMAT OPENPMX_SFORMAT OPENPMX_FFORMAT OPENPMX_SFORMAT "\n", name, lower, value, upper);

			++j;
		}
	}
	
	forcount(i, popmodel.nsigma) {
		if (popmodel.sigmafixed[i] == 0) {
			decode_popmodel_limit(&encode, Hinv, j, z);
			snprintf(upper, sizeof(upper), OPENPMX_FFORMAT, encode.popmodel.sigma[i]);

			decode_popmodel_limit(&encode, Hinv, j, -z); 
			snprintf(lower, sizeof(lower), OPENPMX_FFORMAT, encode.popmodel.sigma[i]);

			let value = popmodel.sigma[i];
			let ind = (source->data._offset1) ? (i + 1) : (i); 
			snprintf(name, sizeof(name), "SIGMA(%i)", ind); 
			info(outstream, OPENPMX_SFORMAT OPENPMX_SFORMAT OPENPMX_FFORMAT OPENPMX_SFORMAT "\n", name, lower, value, upper);

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
					decode_popmodel_limit(&encode, Hinv, j, z);
					snprintf(upper, sizeof(upper), OPENPMX_FFORMAT, encode.popmodel.omega[r][c]);

					decode_popmodel_limit(&encode, Hinv, j, -z); 
					snprintf(lower, sizeof(lower), OPENPMX_FFORMAT, encode.popmodel.omega[r][c]);

					let value = popmodel.omega[r][c];
					let ind = (source->data._offset1) ? (r + 1) : (r); 
					snprintf(name, sizeof(name), "OMEGA(%i,%i)", ind, ind); 
					info(outstream, OPENPMX_SFORMAT OPENPMX_SFORMAT OPENPMX_FFORMAT OPENPMX_SFORMAT "\n", name, lower, value, upper);
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

	/* Eigendecomposition of the correlations in the approximate
	 * covariance matrix, the Hessian inverse. We can skip the 4 (or 2) 
	 * multiplier because we normalize the diagonal to 1 to obtain the
	 * correlation matrix */
	var Heigen = gsl_matrix_alloc(nparam, nparam);
	gsl_matrix_memcpy(Heigen, Hinv);
	scale_to_match_diagonal(Heigen, 0); /* scales for 1 on diagonal */
	var eval = gsl_vector_alloc(nparam);
	var evec = gsl_matrix_alloc(nparam, nparam);
	var w    = gsl_eigen_symmv_alloc(nparam);
	gsl_eigen_symmv(Heigen, eval, evec, w);
	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

	info(outstream, "cov eigendecomposition of correlations in approximate Hessian (J'J) inverse\n");
	info(outstream, OPENPMX_SFORMAT OPENPMX_SFORMAT "\n", "INDEX", "EIGENVALUE");
	forcount(i, nparam)
		info(outstream, OPENPMX_IFORMAT OPENPMX_FFORMAT "\n", i + 1, gsl_vector_get(eval, i));

	let lambda_min = gsl_vector_get(eval, 0);
	let lambda_max = gsl_vector_get(eval, nparam - 1);
	if (lambda_min != 0.)
		info(outstream, "cov condition number %g\n", fabs(lambda_max / lambda_min));

	gsl_eigen_symmv_free(w);
	gsl_matrix_free(evec);
	gsl_vector_free(eval);
	gsl_matrix_free(Heigen);

	/* fallthrough to cleanup */

	assert(j == encode.nparam);

failed:
	if (H)
		gsl_matrix_free(H);
    if (Hinv)
		gsl_matrix_free(Hinv);
    if (outstream)
		fclose(outstream);
		
	free(warm_etas);
	pmx_release_state(&pmx);
}
