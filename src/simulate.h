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

#ifndef OPENPMX_SIMULATE_H
#define OPENPMX_SIMULATE_H

#include "openpmx.h"
#include "idata.h"
#include "popmodel.h"
#include "nonzero.h"
#include "options.h"
#include "scatter.h"

#include <gsl/gsl_rng.h>

#ifdef __cplusplus
extern "C" {
#endif

/* --------------------------------------------------------------------*/

void simulate_with_error_thread(INDIVID* const individ,
								const ADVANFUNCS* const advanfuncs,
								const POPMODEL* const popmodel,
								const NONZERO* const nonzero,
								const OPTIONS* const options,
								const SCATTEROPTIONS* const scatteroptions);

void idata_resample_err(const IDATA* const idata,
						const POPMODEL* const popmodel,
						gsl_rng * const rng);
void idata_resample_eta(IDATA* const idata,
						const POPMODEL* const popmodel,
						gsl_rng * const rng);

#ifdef __cplusplus
}
#endif

#endif
