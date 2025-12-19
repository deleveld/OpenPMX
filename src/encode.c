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
 
/// This file performs the encoding and decoding of a POPMODEL to and 
/// from an unbounded vector. This is a part of the Stage 2 of 
/// optimization.

#include <assert.h>
#include <math.h>
#include <float.h>

#include "encode.h"
#include "linalg.h"
#include "popmodel.h"
#include "utils/c22.h"

#include <gsl/gsl_blas.h>

//#define ESTIMATE_SD_NOT_VAR

//#define ENCODE_SPECIAL
#define ENCODE_ATANH
//#define ENCODE_LINEAR
//#define ENCODE_NONMEM

static int encode_nparam(const POPMODEL* const popmodel,
						 const OMEGAINFO* const omegainfo)
{
	var nparam = 0;

	let ntheta = popmodel->ntheta;
	let thetaestim = popmodel->thetaestim;
	forcount(i, ntheta) {
		if (thetaestim[i] != FIXED)
			++nparam;
	}
	let nsigma = popmodel->nsigma;
	let sigmafixed = popmodel->sigmafixed;
	forcount(i, nsigma) {
		if (sigmafixed[i] == 0)
			++nparam;
	}
	let nonfixed = &omegainfo->nonfixed;
	var ndim = nonfixed->n;
	if (ndim) {
		var rowcol = nonfixed->rowcol;
		forcount(i, ndim) {
			forcount(j, i+1) {
				let r = rowcol[i];
				let c = rowcol[j];
				let fixed = popmodel->omegafixed[r][c];
				if (fixed == 0)
					++nparam;
			}
		}
	}
	return nparam;
}

ENCODE encode_init(const POPMODEL* const popmodel)
{
	var temp_popmodel = *popmodel;
	temp_popmodel.result = (PMXRESULT) { .objfn = DBL_MAX,
										 .type = OBJFN_INVALID,
										 .nparam = 0,
										 .neval = 0 };	
	var temp_omegainfo = omegainfo_init(popmodel->nomega, popmodel->omega, popmodel->omegafixed);
	let n = encode_nparam(popmodel, &temp_omegainfo);

	return (ENCODE) {
		.popmodel = temp_popmodel,
		.omegainfo = temp_omegainfo,
		.nparam = n,
		.offset = { },
		.has_offsets = false,
	};
}

#ifdef ENCODE_SPECIAL
static inline double backward(const double x)
{
	double v = x / 2.;
	double f = v;
	if (v > 0.9)
		f = 0.9 + 0.1*tanh(10.*v - 9.);
	if (v < -0.9)
		f = -0.9 - 0.1*tanh(-10.*v - 9.);
	return f;
}

static inline double forward(const double v)
{
	double f = v;
	if (v > 0.9)
		f = 0.9 - 0.1 * atanh(9. - 10.*v);
	if (v < -0.9)
		f = -0.9 + 0.1 * atanh(10.*v + 9.);
	return f * 2.;
}
#endif

static double theta_transform(const double value, const double lower, const double upper)
{
#ifdef ENCODE_SPECIAL
	return forward(2. * (value - lower) / (upper - lower) - 1.);
#endif
#ifdef ENCODE_ATANH
	return atanh(2. * (value - lower) / (upper - lower) - 1.);
#endif
#ifdef ENCODE_LINEAR
	if (value < lower)
		return -1.;
	if (value > lower)
		return 1.;
	return 2. * (value - lower) / (upper - lower) - 1.;
#endif
#ifdef ENCODE_NONMEM
// Same as NONMEM
	return log(value - lower) - log(upper - value);
#endif
}

static double theta_untransform(const double f, const double lower, const double upper)
{
#ifdef ENCODE_SPECIAL
	return lower + (backward(f) + 1.) * (upper - lower) / 2.;
#endif
#ifdef ENCODE_ATANH
	return lower + (tanh(f) + 1.) * (upper - lower) / 2.;
#endif
#ifdef ENCODE_LINEAR
	return lower + (f + 1.) * (upper - lower) / 2.;
#endif
#ifdef ENCODE_NONMEM
// Same as NONMEM
// https://www.wolframalpha.com/input?i2d=true&i=solve+for+x%3B+w+%3D+log%5C%2840%29x-l%5C%2841%29-log%5C%2840%29u-x%5C%2841%29
	return (lower + upper * exp(f)) / (exp(f) + 1.);
#endif
}

/* scale the matrix to have the same diagonal as the reference
 * If ref is NULL then assume identity matrix */
static void scale_to_match_diagonal(gsl_matrix* matrix, const gsl_matrix* ref)
{
	let n = matrix->size1;

	double tempdata[OPENPMX_OMEGA_MAX * OPENPMX_OMEGA_MAX];
	var temp = gsl_matrix_view_array(tempdata, n, n);

	double scaledata[OPENPMX_OMEGA_MAX * OPENPMX_OMEGA_MAX];
	var scale = gsl_matrix_view_array(scaledata, n, n);
	gsl_matrix_set_zero(&scale.matrix);

	forcount(i, n) {
		var s = 1.;
		if (ref)
			s = gsl_matrix_get(ref, i, i);
		var v = gsl_matrix_get(matrix, i, i);
		double x;
		if (v == 0)
			x = 0.;
		else
			x = sqrt(s/v);
		gsl_matrix_set(&scale.matrix, i, i, x);
	}
	/* correct scale as S*omega*ST */
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., &scale.matrix, matrix, 0., &temp.matrix);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., &temp.matrix, &scale.matrix, 0., matrix);
}

/// Calling encode_offset() makes such that a zero-initialized vector 
/// will reproduce the given POPMODEL. 

/* encode the popmodel, setting the offsets */
void encode_offset(ENCODE* const encode, const POPMODEL* const popmodel)
{
	encode->popmodel = *popmodel;

	/* foce advan need to see updated omega cholesky and lndet */
	/* we assume the empty rows and cols dont change */
	var omegainfo = &encode->omegainfo;
	omegainfo_update_inverse_lndet(omegainfo, popmodel->omega);

	/* we write directly into the offset */
	var x = encode->offset;
	encode->has_offsets = true;
	
/// Different transformations are possible to address theta bounds but 
/// the default is to use tanh/atanh. 
	var n = 0;
	let etheta = popmodel->theta;
	let ntheta = popmodel->ntheta;
	let thetaestim = popmodel->thetaestim;
	forcount(i, ntheta) {
		if (thetaestim[i] != FIXED) {
			let l = popmodel->lower[i];
			let v = etheta[i];
			let u = popmodel->upper[i];
			x[n] = theta_transform(v, l, u);
			++n;
		}
	}

/// Sigma is log-transformed. 
	let sigma = popmodel->sigma;
	let nsigma = popmodel->nsigma;
	let sigmafixed = popmodel->sigmafixed;
	forcount(i, nsigma) {
		if (sigmafixed[i] == 0) {
#ifdef ESTIMATE_SD_NOT_VAR
			let s = sqrt(sigma[i]);
#else
			let s = sigma[i];
#endif
			x[n] = log(s);
			++n;
		}
	}

/// The diagonal of omega is log-transformed.
///
/// The off-diagonals of omega are encoded as a specialized Cholesky decomposition of a 
/// correlation matrix. Only the lower triangular is encoded. The 
/// approach is described in: <https://mc-stan.org/docs/reference-manual/transforms.html#cholesky-factors-of-correlation-matrices>
/// and also: <https://mc-stan.org/docs/reference-manual/transforms.html#correlation-matrix-transform.section>
	let nonfixed = &omegainfo->nonfixed;
	var ndim = nonfixed->n;
	if (ndim) {
		var cholesky = gsl_matrix_alloc(ndim, ndim);

		/* make a correlation matrix and take the cholesky of this */
		/* because we scale matrix to 1 on diagonal it becomes a correlation matrix */
		reduced_omega_init(cholesky, popmodel->omega, nonfixed->rowcol, nonfixed->n);
		scale_to_match_diagonal(cholesky, 0);
		cholesky_decomposition(cholesky);
		let fail = gsl_matrix_transpose(cholesky);
		assert(!fail);

		/* https://mc-stan.org/docs/reference-manual/transforms.html#cholesky-factors-of-correlation-matrices */
		/* The next step from the Cholesky factor w back to the array z of canonical partial correlations
		 * (CPCs) is simplified by the ordering of the elements in the definition of w, which when inverted yields" */
		var z = gsl_matrix_alloc(ndim, ndim);
		gsl_matrix_set_all(z, 0.);
		forcount(j, ndim) {
			forcount(i, ndim) {
				var v = 0.;
				if (i == 0 && i < j) {
					v = gsl_matrix_get(cholesky, i, j);
				} else if (i > 0 && i < j) {
					v = gsl_matrix_get(cholesky, i, j);
					for (int ip=0; ip<i; ip++) {
						let zipj = gsl_matrix_get(cholesky, ip, j);
						v *= pow(1. - zipj*zipj, -0.5);
					}
				}
				gsl_matrix_set(z, i, j, v);
			}
		}
		var rowcol = nonfixed->rowcol;
		forcount(i, ndim) {
			forcount(j, i+1) {
				let r = rowcol[i];
				let c = rowcol[j];
				let fixed = popmodel->omegafixed[r][c];
				if (fixed == 0) {
					if (i == j) {
#ifdef ESTIMATE_SD_NOT_VAR
						let s = sqrt(popmodel->omega[r][c]);
#else
						let s = popmodel->omega[r][c];
#endif
						x[n] = log(s);
					} else {
						let s = gsl_matrix_get(z, j, i);
						/* see above webpage close to the text
						 * "The final stage of the transform reverses the
						 * hyperbolic tangent transform, which is defined by" */
						x[n] = atanh(s);
					}
					++n;
				}
			}
		}
		gsl_matrix_free(z);
		gsl_matrix_free(cholesky);
	}
	assert(n == encode->nparam);
}

/// If an OMEGA_SAME block is used then (OMEGASAME() in openpmxtran)
/// then omega is filled by block transposed by its size. Basically it
/// reproduces a block of omega by a lower block. This allows two eta
/// values with the same variances to be estimated.
static void fill_in_OMEGA_SAME_blocks(POPMODEL* const popmodel)
{
	/* correct for OMEGA_SAME blocks */
	let nblock = popmodel->nblock;
	var d = 0;
	forcount(i, nblock) {
		let btype = popmodel->blocktype[i];
		let n = popmodel->blockdim[i];
		if (btype == OMEGA_SAME) {
			forcount(r, n) {
				for (var c=0; c<r; c++) {
					let v = popmodel->omega[d + c - n][d + r - n];
					popmodel->omega[d + c][d + r] = v;
					popmodel->omega[d + r][d + c] = v;
				}
				let v = popmodel->omega[d + r - n][d + r - n];
				popmodel->omega[d + r][d + r] = v;
			}
		}
		d += n;
	}
}

static void popmodel_omega_update(POPMODEL* const popmodel,
								  const gsl_matrix * const omegareduced_nonfixed,
								  const int * const rowcol)
{
	/* have to expand the reduced omega to fill in the original omega */
	let n = omegareduced_nonfixed->size1;
	forcount(i, n) {
		forcount(j, n) {
			let r = rowcol[i];
			let c = rowcol[j];
			var v = gsl_matrix_get(omegareduced_nonfixed, i, j);
			popmodel->omega[r][c] = v;
		}
	}
	fill_in_OMEGA_SAME_blocks(popmodel);
}

void encode_update(ENCODE* encode, const double* x)
{
	assert(encode->has_offsets == true);

	var popmodel = &encode->popmodel;
	var omegainfo = &encode->omegainfo;
	let offset = encode->offset;
	
	var n = 0;
	let etheta = popmodel->theta;
	let ntheta = popmodel->ntheta;
	let thetaestim = popmodel->thetaestim;
	forcount(i, ntheta) {
		if (thetaestim[i] != FIXED) {
			var l = popmodel->lower[i];
			var v = (x[n] + offset[n]);
			var u = popmodel->upper[i];
			var estvalue = theta_untransform(v, l, u);
			etheta[i] = estvalue;
			++n;
		}
	}
	let sigma = popmodel->sigma;
	let nsigma = popmodel->nsigma;
	let sigmafixed = popmodel->sigmafixed;
	forcount(i, nsigma) {
		if (sigmafixed[i] == 0) {
			let s = (x[n] + offset[n]);
#ifdef ESTIMATE_SD_NOT_VAR
			sigma[i] = pow(exp(s), 2);
#else
			sigma[i] = exp(s);
#endif
			++n;
		}
	}

	let nonfixed = &omegainfo->nonfixed;
	var ndim = nonfixed->n;
	if (ndim) {
		var rowcol = nonfixed->rowcol;
		var S = gsl_matrix_alloc(ndim, ndim);
		var z = gsl_matrix_alloc(ndim, ndim);
		gsl_matrix_set_zero(S);
		gsl_matrix_set_zero(z);

		/* copy over the diagonal so we have any fixed variances set */
		forcount(i, ndim) {
			let rc = rowcol[i];
			let v = popmodel->omega[rc][rc];
			gsl_matrix_set(S, i, i, v);
		}

		forcount(i, ndim) {
			forcount(j, i+1) {
				/* parameters we estimated we extract from the vector */
				let r = rowcol[i];
				let c = rowcol[j];
				let fixed = popmodel->omegafixed[r][c];
				if (fixed == 0) {
					let s = (x[n] + offset[n]);
					if (i == j) {
#ifdef ESTIMATE_SD_NOT_VAR
						gsl_matrix_set(S, i, i, pow(exp(s), 2));
#else
						gsl_matrix_set(S, i, i, exp(s));
#endif
					} else
						gsl_matrix_set(z, j, i, tanh(s));
					++n;
				}
			}
		}

		/* https://mc-stan.org/docs/reference-manual/transforms.html#correlation-matrix-transform.section */
		/* near the text "In Stan, the LKJ transform is reformulated in terms of a Cholesky
		 * factor w of the final correlation matrix, defined for 1≤i,j≤K by */
		var cholesky = gsl_matrix_alloc(ndim, ndim);
		gsl_matrix_set_zero(cholesky);
		forcount(i, ndim) {
			forcount(j, ndim) {
				var wij = 0.;
				if (i == 0 && j == 0) {
					wij = 1.;
				} else if (i > 0 && i == j) {
					var p = 1.;
					for (int ip=0; ip<i; ip++) {
						let zipj = gsl_matrix_get(z, ip, j);
						p *= pow(1. - zipj*zipj, 0.5);
					}
					wij = p;
				} else if (i == 0 && i < j) {
					wij = gsl_matrix_get(z, i, j);
				} else if (i > 0 && i < j) {
					let zij = gsl_matrix_get(z, i, j);
					double p = zij;
					for (int ip=0; ip<i; ip++) {
						let zipj = gsl_matrix_get(z, ip, j);
						p *= pow(1. - zipj*zipj, 0.5);
					}
					wij = p;
				}
				gsl_matrix_set(cholesky, i, j, wij);
			}
		}
		var temp = gsl_matrix_alloc(ndim, ndim);
		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., cholesky, cholesky, 0., temp);

		scale_to_match_diagonal(temp, S);
		gsl_matrix_free(S);
		gsl_matrix_free(cholesky);
		gsl_matrix_free(z);

		/* update full omega in the popmodel. We have to update the
		 * OMEGAINFO later so that the objective function log(det()) term
		 * is correct */
		popmodel_omega_update(popmodel, temp, omegainfo->nonfixed.rowcol);
		gsl_matrix_free(temp);
	}
	assert(n == encode->nparam);

	/* foce advan need to see updated omega cholesky and lndet */
	/* we assume the empty rows and cols dont change */
	omegainfo_update_inverse_lndet(omegainfo, popmodel->omega);

	/* invalidate the objfn because we dont know it anymore */
	popmodel->result = (PMXRESULT) { .objfn = DBL_MAX,
									 .type = OBJFN_INVALID,
									 .nparam = 0,
									 .neval = 0 };
}

