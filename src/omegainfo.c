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

#include "omegainfo.h"
#include "linalg.h"
#include "utils/c22.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "openpmx_compile_options.h"

void reduced_omega_init(gsl_matrix* matrix,
						const double omega[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX],
						const int* const rowcol,
						const int n)
{
	assert(n == (int)matrix->size1);
	forcount(i, n) {
		forcount(j, n) {
			let r = rowcol[i];
			let c = rowcol[j];
			let v = omega[r][c];
			gsl_matrix_set(matrix, i, j, v);
		}
	}
}

void omegainfo_update_inverse_lndet(OMEGAINFO* const omegainfo,
									const double omega[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX])
{
	assert(omegainfo->nonzero.n >= 0);
	let n = omegainfo->nonzero.n;
	var cholesky = gsl_matrix_view_array(omegainfo->nonzero.choleskydata, n, n);
	var inverse = gsl_matrix_view_array(omegainfo->nonzero.inversedata, n, n);

	/* The reduced non-zero matrix is used during stage 1 estimation.
	 * These are the etas that we estimate in stage 1 and it log(det)
	 * is a term in the individual objective function. */
	if (n) {
		/* get the nonzero omega matrix to make the cholesky which we
		 * use to find inverse and log(det()) */
		reduced_omega_init(&cholesky.matrix, omega, omegainfo->nonzero.rowcol, omegainfo->nonzero.n);
		cholesky_decomposition(&cholesky.matrix, "nonzero");

		/* fill in data of nonzero inverse matrix */
		gsl_matrix_memcpy(&inverse.matrix, &cholesky.matrix);
		gsl_linalg_cholesky_invert(&inverse.matrix);

		/* log(det(omega)) term is a part of the objective function */
		let lndet = matrix_lndet_from_cholesky(&cholesky.matrix);
		assert(isfinite(lndet) == 1);
		omegainfo->omega_nonzero_lndet = lndet;

	} else {
		omegainfo->omega_nonzero_lndet = 0.;
		gsl_matrix_set_zero(&cholesky.matrix);
		gsl_matrix_set_zero(&inverse.matrix);
	}
}

OMEGAINFO omegainfo_init(const int nomega,
						 const double omega[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX],
						 const int omegafixed[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX])
{
	OMEGAINFO ret = { 0 };
	ret.omega_nonzero_lndet = 0.;
	ret.nonzero.n = 0;
	ret.nonfixed.n = 0;

	var maxdiag = -DBL_MAX;
	forcount(i, nomega) {
		let v = omega[i][i];
		if (v > maxdiag)
			maxdiag = v;
	}

	if (maxdiag > 0.) {
		forcount(i, nomega) {

			/* sum of this row and coloum */
			var sum = 0.;
			forcount(j, i+1) {
				sum += omega[i][j];
				sum += omega[j][i];
			}

			/* there is something non-zero in this row/col */
			if (sum != 0.) {
				ret.nonzero.rowcol[ret.nonzero.n] = i;
				++ret.nonzero.n;
			}

			/* if the diagonal is SAME then the whole row/col is not going
			 * to be estimated. the negative variances (FIXED) on the
			 * diagonal are handeled during decode */
			let fixed = omegafixed[i][i];
			if (sum != 0. && fixed != 2) {
				ret.nonfixed.rowcol[ret.nonfixed.n] = i;
				++ret.nonfixed.n;
			}
		}
	}

	/* fill in the nonzero and nonfixed matricies with the right rows and
	 * coloumns from omega and then fill in the cholesky decomposition */
	omegainfo_update_inverse_lndet(&ret, omega);

	return ret;
}

