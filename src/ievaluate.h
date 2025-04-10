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

#ifndef OPENPMX_EVALUATE_H
#define OPENPMX_EVALUATE_H

#include <stdlib.h>
#include <math.h>

#include "openpmx.h"

#ifdef __cplusplus
extern "C" {
#endif

/* --------------------------------------------------------------------*/

typedef struct {
	const RECORD* const record;
	const int nrecord;
	const ADVANFUNCS* const advanfuncs;
	const POPPARAM popparam;
	FILE* logstream;
} IEVALUATE_ARGS;

double individual_fasteval(const IEVALUATE_ARGS* const ievaluate_args);
void individual_evaluate(const IEVALUATE_ARGS* const ievaluate_args,
						 IMODEL* const imodel_saved,
						 PREDICTVARS* const predictvars_saved,
						 double* const istate,
						 double* const YHAT,
						 double* const YHATVAR,
						 double* const obs_lndet,
						 double* const obs_min2ll);
void individual_checkout(const IEVALUATE_ARGS* const ievaluate_args);
void individual_simulate(const IEVALUATE_ARGS* const ievaluate_args,
						 IMODEL* const imodel_saved,
						 PREDICTVARS* const predictvars_saved,
						 double* const istate,
						 double* const individ_yhat,
						 double* const individ_yhatvar,
						 double* const individ_pred,
						 const double* const isimerr);

#ifdef __cplusplus
}
#endif

#endif
