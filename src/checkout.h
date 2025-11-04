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

#ifndef OPENPMX_CHECKOUT_H
#define OPENPMX_CHECKOUT_H

#include <stdio.h>

#include "idata.h"
#include "options.h"
#include "openpmx.h"

#ifdef __cplusplus
extern "C" {
#endif

void idata_checkout(const IDATA* const idata,
					const ADVANFUNCS* const advanfuncs,
					const POPMODEL* const popmodel,
					const OPTIONS* const options,
					FILE* logstream);

#ifdef __cplusplus
}
#endif

#endif
