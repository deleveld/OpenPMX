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

#ifndef PMXSTATE_H
#define PMXSTATE_H

#include "popmodel.h"
#include "idata.h"

#include <gsl/gsl_rng.h>

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------*/

typedef struct PMXSTATE {
	const ADVANFUNCS* const advanfuncs;
	IDATA idata;
	int tablecount;
	gsl_rng* rng;
} PMXSTATE;

void pmx_update_from_popmodel(OPENPMX* const pmx, const POPMODEL* const popmodel);

void pmxstate_ensure(OPENPMX* const pmx);

#ifdef __cplusplus
}
#endif

#endif
