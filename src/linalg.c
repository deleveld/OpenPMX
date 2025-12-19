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
 
/// This file contains linear algebra functions for Cholesky 
/// decomposition. This is necessary for calculating log(det()) of a 
/// covariance matrix, and calculating the sample likelihood.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "linalg.h"
#include "utils/c22.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

int cholesky_decomposition(gsl_matrix* matrix)
{
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

double sample_min2ll_from_cholesky(const double* const data,
								   const gsl_matrix* const chol)
{
	let n = chol->size1;
	let x = gsl_vector_const_view_array(data, n);

	/* Basically this 'decorrelates' the sample with respect to the covariance matrix
	 (via its cholesky decomposition) which gives the independent samples. The sum of
	 the square is the likelihood */
	double rdata[n];
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

double sample_min2ll_from_inverse(const double* const data,
								  const gsl_matrix* const inverse)
{
	let n = inverse->size1;
	double ydata[n];
	var y = gsl_vector_view_array(ydata, n);
	let x = gsl_vector_const_view_array(data, n);

	gsl_blas_dgemv(CblasNoTrans, 1., inverse, &x.vector, 0., &y.vector);

//	var lik = 0.;
//	gsl_blas_ddot(&x.vector, &y.vector, &lik);
//	return lik;
	double sum = 0.;
	forcount(i, n) {
		let _x = gsl_vector_get(&x.vector, i);
		let _y = gsl_vector_get(&y.vector, i);
		sum += _x * _y;
	}
	return sum;
}


