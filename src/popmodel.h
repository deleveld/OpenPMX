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

#ifndef OPENPMX_POPMODEL_H
#define OPENPMX_POPMODEL_H

#include <stdio.h>

#include "openpmx.h"
#include "utils/errctx.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef typeof(((OPENPMX){0}).result) PMXRESULT;

typedef struct {
	double lower[OPENPMX_THETA_MAX];
	double theta[OPENPMX_THETA_MAX];
	double upper[OPENPMX_THETA_MAX];
	unsigned char thetaestim[OPENPMX_THETA_MAX]; 		/* openpmx.theta[0].type */
	int ntheta;

	unsigned char blocktype[OPENPMX_OMEGABLOCK_MAX];	/* openpmx.omega[0].type */
	int blockdim[OPENPMX_OMEGABLOCK_MAX];
	int nblock;

	double omega[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX];
	unsigned char omegafixed[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX];
	int nomega;

	double sigma[OPENPMX_SIGMA_MAX];
	unsigned char sigmafixed[OPENPMX_SIGMA_MAX];
	int nsigma;

	/* filled in during estimation iterations */
	PMXRESULT result;
} POPMODEL;

typedef typeof(((OPENPMX){0}).theta[0]) THETATYPE;
typedef typeof(((OPENPMX){0}).omega[0]) OMEGABLOCKSTYPE;

POPMODEL popmodel_init(const THETATYPE* const theta,
					   const OMEGABLOCKSTYPE* const omegablocks,
					   const double* sigma,
					   ERRCTX* errctx);

void extfile_header(FILE * f,
					const POPMODEL* const popmodel,
					const bool _offset1);
void extfile_append(FILE* f,
					const POPMODEL* const popmodel,
					const double runtime_s,
					const int iter,
					const int ineval);
void extfile_trailer(FILE* f, const POPMODEL* const popmodel,
					 const double runtime_s,
					 const int ineval);

void popmodel_eval_information(const POPMODEL* const popmodel,
							   const double runtime_s,
							   const int neval,
							   const bool details,
							   FILE* outstream,
							   FILE* extstream);

void popmodel_information(FILE* f2,
						  const POPMODEL* const popmodel,
						  const double timestamp);


#ifdef __cplusplus
}
#endif

#endif
