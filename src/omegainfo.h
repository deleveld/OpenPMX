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

#ifndef OPENPMX_OMEGAINFO_H
#define OPENPMX_OMEGAINFO_H

#include "openpmx.h"
#include "nonzero.h"

#include <gsl/gsl_matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	NONZERO nonzero;
	struct {
		int rowcol[OPENPMX_OMEGA_MAX];
		int n;
	} nonfixed;
	double omega_nonzero_lndet;
} OMEGAINFO;

OMEGAINFO omegainfo_init(const int nomega,
						 const double omega[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX],
						 const int omegafixed[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX]);

void omegainfo_update_inverse_lndet(OMEGAINFO* const omegainfo,
									const double omega[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX]);

void reduced_omega_init(gsl_matrix* matrix,
						const double omega[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX],
						const int* const rowcol,
						const int n);
 
#ifdef __cplusplus
}
#endif

#endif
