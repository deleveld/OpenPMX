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

#ifndef OPENPMX_EVALUATE_H
#define OPENPMX_EVALUATE_H

#include <stdlib.h>
#include <math.h>

#include "openpmx.h"

#ifdef __cplusplus
extern "C" {
#endif

/* --------------------------------------------------------------------*/
#define USE_KAHAN_SUM

/* use Kahan summation to reduce numerical issues in summing the
 * objective function calculation some methods use gradient from
 * objective function to this needs to be as accurate as possible
 * in the context of rounding errors.  It is not really clear whether
 * this is affected by compiler optimization. This is the reason why
 * the KAHAN pointer is volatile */
/* https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements */
typedef struct {
#ifdef USE_KAHAN_SUM
	double _c;
	double _sum;
#else
	double _sum;
#endif
} KAHAN;

static inline void KAHAN_ADD(const double val, volatile KAHAN* const _k)
{
#ifdef USE_KAHAN_SUM
	const double t = _k->_sum + val;
	if (fabs(_k->_sum) >= fabs(val))
		_k->_c += (_k->_sum - t) + val;
	else
		_k->_c += (val - t) + _k->_sum;
	_k->_sum = t;
#else
	_k->_sum += val;
#endif
}

static inline double KAHAN_SUM(const KAHAN* const _k)
{
#ifdef USE_KAHAN_SUM
	return _k->_sum + _k->_c;
#else
	return _k->_sum;
#endif
}

typedef struct {
	const RECORD* const record;
	const int nrecord;
	const ADVANFUNCS* const advanfuncs;
	const POPPARAM popparam;
	FILE* logstream;
} IEVALUATE_ARGS;

void individual_evaluate(const IEVALUATE_ARGS* const ievaluate_args,
						 IMODEL* const imodel_saved,
						 PREDICTVARS* const predictvars_saved,
						 double* const istate,
						 double* const YHAT,
						 double* const YHATVAR,
						 KAHAN* const obs_lndet,
						 KAHAN* const obs_min2ll);
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
