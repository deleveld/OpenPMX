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

#include "nonzero.h"
#include "linalg.h"

/* different ways to calculate sample min2ll */
//#define SAMPLE_MIN2LL_FROM_INVERSE
#define SAMPLE_MIN2LL_FROM_CHOLESKY

void reduce_eta(double* reta,
				const double fulleta[static OPENPMX_OMEGA_MAX],
				const NONZERO* const nonzero)
{
	for (int i=0; i<nonzero->n; i++) {
		const int col = nonzero->rowcol[i];
		reta[i] = fulleta[col];
	}
}

void unreduce_eta(double fulleta[static OPENPMX_OMEGA_MAX],
				  const double* reta,
				  const NONZERO* const nonzero)
{
	memset(fulleta, 0, OPENPMX_OMEGA_MAX * sizeof(double));
	for (int i=0; i<nonzero->n; i++) {
		const int col = nonzero->rowcol[i];
		const double v = reta[i];
		fulleta[col] = v;
	}
}

/* functions for Bae and Yim objective function Term 3 */
double sample_min2ll(const int nreta,
					 const double reta[static nreta],
					 const NONZERO* const nonzero)
{
	assert(nreta == nonzero->n);
	double term3 = 0.;
#ifdef SAMPLE_MIN2LL_FROM_INVERSE
	const double* omegainverse_data = nonzero->inversedata;
	const gsl_matrix_const_view  = gsl_matrix_const_view_array(omegainverse_data, nreta, nreta);
	term3 += sample_min2ll_from_inverse(reta, &omegainverse.matrix);
#endif
#ifdef SAMPLE_MIN2LL_FROM_CHOLESKY
	const double* omegacholesky_data = nonzero->choleskydata;
	const gsl_matrix_const_view omegacholesky = gsl_matrix_const_view_array(omegacholesky_data, nreta, nreta);
	term3 += sample_min2ll_from_cholesky(reta, &omegacholesky.matrix);
#endif

#ifdef SAMPLE_MIN2LL_FROM_INVERSE
#ifdef SAMPLE_MIN2LL_FROM_CHOLESKY
	term3 /= 2.;
#endif
#endif
//	assert(gsl_finite(term3) == 1);
	return term3;
}

