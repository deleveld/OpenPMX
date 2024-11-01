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

#ifndef OPENPMX_STAGE1_H
#define OPENPMX_STAGE1_H

#include "scatter.h"
#include "popmodel.h"

#ifdef __cplusplus
extern "C" {
#endif

void stage1_thread(INDIVID* const individ,
				   const ADVANFUNCS* const advanfuncs,
				   const POPMODEL* const popmodel,
				   const NONZERO* const nonzero,
				   const OPTIONS* const options,
				   const SCATTEROPTIONS* const scatteroptions);

#ifdef __cplusplus
}
#endif

#endif
