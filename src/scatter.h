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

#ifndef OPENPMX_SCATTER_H
#define OPENPMX_SCATTER_H

#include <stdio.h>

#include "idata.h"
#include "popmodel.h"
#include "nonzero.h"
#include "openpmx_internal.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	bool stage1_order;
	bool checkout_errors;
	FILE* logstream;
} SCATTEROPTIONS;

typedef void (*THREADTASK)(INDIVID* const individ,
						   const ADVANFUNCS* const advanfuncs,
						   const POPMODEL* const popmodel,
						   const NONZERO* const nonzero,
						   const OPTIONS* const options,
						   const SCATTEROPTIONS* const scatteroptions);
						   
void scatter_threads(const IDATA* const idata,
					 const ADVANFUNCS* const advanfuncs,
					 const POPMODEL* const popmodel,
					 const NONZERO* const nonzero,
					 const OPTIONS* const options,
					 SCATTEROPTIONS* scatteroptions,
					 THREADTASK threadtask);

void scatter_cleanup(const bool progress);

#ifdef __cplusplus
}
#endif

#endif
