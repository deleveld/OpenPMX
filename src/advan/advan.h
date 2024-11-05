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

#ifndef OPENPMX_ADVAN_INTERNALH
#define OPENPMX_ADVAN_INTERNALH

#include <stdlib.h>
#include <stdio.h>

#include "dataconfig/recordinfo.h"
#include "utils/vector.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ADVAN ADVAN;

typedef struct ADVANFUNCS {
	const int advan_size;
	void (*construct)(ADVAN* advan, const ADVANFUNCS* const advanfuncs);
	void (*destruct)(ADVAN * advan);
	void (*info)(const ADVANFUNCS* const advanfuncs, FILE* f);

	void (*reset)(ADVAN * advan, const int full);
	void (*interval)(ADVAN* advan,
					 const IMODEL* const imodel,
					 const RECORD* const record,
					 double* state,
					 const POPPARAM* const popparam,
					 const double endtime,
					 const double* rates);

	const ADVANCONFIG* const advanconfig;
	const RECORDINFO recordinfo;
	const int nstate;
} ADVANFUNCS;

typedef struct {
	int cmt;
	double amt;
	double rate;
	double start;
	double end;
} ADVANINFUSION;

typedef struct ADVAN {
	const ADVANFUNCS* advanfuncs;

	double time;
	double state[OPENPMX_STATE_MAX];
	double amtlag[OPENPMX_STATE_MAX];
	double bioavail[OPENPMX_STATE_MAX];

	ADVANINFUSION _infusions_buffer[OPENPMX_SIMULINFUSION_MAX];
	VECTOR(ADVANINFUSION) infusions;

	int initcount;
} ADVAN;

void advan_base_construct(ADVAN* advanbase, const ADVANFUNCS* advanfuncs);
void advan_base_destruct(ADVAN* advanbase);

PREDICTSTATE advan_advance(ADVAN* const advan,
						   IMODEL* const imodel,
						   const RECORD* const record,
						   const POPPARAM* const popparam);

ADVANFUNCS* advanfuncs_alloc(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);
void advanfuncs_free(ADVANFUNCS* advan);

/* Used by some ODE methods to send arguments to the user DIFFEQN function */
typedef struct ADVANCER_DIFFEQN_CALLBACK_ARGS {
	ADVAN_DIFFEQN diffeqn;
	const ADVAN* advan;
	const IMODEL* imodel;
	const RECORD* record;
	const POPPARAM* popparam;
	const double* rates;
	int nstate;
} ADVANCER_DIFFEQN_CALLBACK_ARGS;

#ifdef __cplusplus
}
#endif

#endif
