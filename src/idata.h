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

#ifndef OPENPMX_IDATA_H
#define OPENPMX_IDATA_H

#include "dataconfig/dataconfig.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	const double ID;
	const RECORD* const record;	/* view into data in RECORDINFO */
	const int nrecord;
	const int nobs;

	/* the allocated memory is at the first individual */
	IMODEL* const imodel;
	double* const istate;
	double* const eta;
	double* const icov;
	double* const yhat;
	double* const yhatvar;
	double* const pred;
	PREDICTVARS* const predictvars;
	double* icovsample;
	double* icovweight;
	double* isimerr;	/* only for simulation */

	/* objective function */
	double obs_min2ll;
	double obs_lndet;
	double eta_min2ll;
	double icov_lndet;
	double iobjfn;

	/* execution time */
	double eval_msec;
	double stage1_msec;
	int ineval;
} INDIVID;

typedef struct {
	const int ndata;
	const int nindivid;
	const int nobs;
	const int ntheta;
	const int nomega;
	const int nsigma;
	const int nstate;
	const int imodel_size;
	const int predictvars_size;

	INDIVID* const individ;		/* information about each individual */
	/* pointers to these can be found at the first individual
	IMODEL* const imodel;
	double* const istate;
	double* const eta;
	double* const icov;
	double* const yhat;
	double* const yhatvar;
	double* const pred;
	PREDICTVARS* const predictvars;
	double* simerr; */
} IDATA;

IDATA idata_construct(const RECORDINFO* const recordinfo,
					  const int ntheta,
					  const int nomega,
					  const int nsigma,
					  const int nstate,
					  const int imodel_size,
					  const int predictvars_size);
void idata_destruct(IDATA* const idata);

double idata_objfn(const IDATA* const idata,
				   const double omega_nonzero_lndet);
int idata_ineval(const IDATA* const idata, const bool reset);

void idata_reset_eta(IDATA* const idata, const double* eta);

double* idata_alloc_simerr(const IDATA* const idata);
void idata_free_simerr(const IDATA* const idata);

void idata_alloc_icovresample(const IDATA* const idata);
void idata_free_icovresample(const IDATA* const idata);

void table_phi_idata(const char* filename,
					 const IDATA* const idata,
					 const bool _offset1);
void table_icov_resample_idata(const char* filename,
							   const IDATA* const idata,
							   const bool _offset1);

#ifdef __cplusplus
}
#endif

#endif
