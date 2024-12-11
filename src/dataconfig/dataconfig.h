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

#ifndef OPENPMX_DATACONFIG_INTERNAL_H
#define OPENPMX_DATACONFIG_INTERNAL_H

#include <stddef.h>
#include "openpmx.h"

#ifdef __cplusplus
extern "C" {
#endif

int structinfo_find_offset(const char* name, const STRUCTINFO* const structinfo);

#ifdef __cplusplus
}
#endif

#endif
