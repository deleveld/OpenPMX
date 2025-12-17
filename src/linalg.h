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

#ifndef OPENPMX_LINALG_H
#define OPENPMX_LINALG_H

#include <gsl/gsl_matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

int cholesky_decomposition(gsl_matrix* matrix);
double matrix_lndet_from_cholesky(const gsl_matrix* const chol);

double sample_min2ll_from_cholesky(const double* const data,
								   const gsl_matrix* const chol);
double sample_min2ll_from_inverse(const double* const data,
								  const gsl_matrix* const inverse);

 
#ifdef __cplusplus
}
#endif

#endif
