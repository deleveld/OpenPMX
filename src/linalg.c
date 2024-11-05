/* 
 * This file is part of OpenPMX (https://github.com/deleveld/openpmx).
 * Copyright (c) 2022 Douglas Eleveld.
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

#include "linalg.h"
#include "utils/c22.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "openpmx_compile_options.h"

int cholesky_decomposition(gsl_matrix* matrix, const char* matrixname)
{
	/* TODO: finish this */
	(void)matrixname;
	
	let n = matrix->size1;
	assert(n == matrix->size2);

	var mat = gsl_matrix_alloc(n, n);
	gsl_matrix_memcpy(mat, matrix);

//	forcount(i, n) {
//		forcount(j, n) 
//			printf("%f ", gsl_matrix_get(mat, i, j));
//		printf("\n");
//	}

	/* try to do cholesky decomposition */
	let oldhandler = gsl_set_error_handler_off();
	let err = gsl_linalg_cholesky_decomp1(matrix);
	gsl_set_error_handler(oldhandler);

	/* error about not positive definite */
	assert(err != GSL_EDOM);

	/* make the cholesky only lower diagonal */
	forcount(i, n)
		forcount(j, i)
			gsl_matrix_set(matrix, j, i, 0.);

	gsl_matrix_free(mat);

	return 0;
}

double matrix_lndet_from_cholesky(const gsl_matrix* const chol)
{
	/* the det of a matrix is the 2 time the product of the diagonal of the cholesky
	   https://math.stackexchange.com/questions/2001041/logarithm-of-the-determinant-of-a-positive-definite-matrix
	   so we can sum the logs along the diagonal. */
	let n = chol->size1;

	double sum = 0.;
	forcount(i, n) {
		let v = 2. * log(gsl_matrix_get(chol, i, i));
		sum += v;
	}
	return sum;
}

