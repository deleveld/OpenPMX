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

#ifndef OPENPMX_INTERNAL_H
#define OPENPMX_INTERNAL_H

#include "popmodel.h"
#include "idata.h"

#include <gsl/gsl_rng.h>

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------*/

typedef struct {
	int nthread;
	const bool _offset1;
	bool details;
	bool verbose;
	ESTIMCONFIG estimate;
	SIMCONFIG simulate;
} OPTIONS;

OPTIONS options_default(const OPTIONS* const opt1);

STAGE1CONFIG stage1config_default(const STAGE1CONFIG* const stage1);
SIMCONFIG simconfig_default(const SIMCONFIG* const simulate);
ESTIMCONFIG estimconfig_default(const ESTIMCONFIG* const estimate);
OPTIONS options_init_from_pmx(const OPENPMX* const pmx);

POPPARAM popparam_init(const POPMODEL* const popmodel,
					   const ADVANFUNCS* const advanfuncs,
					   const double eta[static OPENPMX_OMEGA_MAX]);

typedef struct PMXSTATE {
	const ADVANFUNCS* const advanfuncs;
	IDATA idata;
	int tablecount;
	gsl_rng* rng;
} PMXSTATE;

void pmxstate_ensure(OPENPMX* const pmx);

#ifdef __cplusplus
}
#endif

#endif
