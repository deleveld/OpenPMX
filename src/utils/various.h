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

#include <stdio.h>
#include <time.h>

#ifndef OPENPMX_UTIL_H
#define OPENPMX_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

FILE* results_fopen(const char* name, const char* ext, const char* mode);

double timespec_time_difference(const struct timespec* const begin,
								const struct timespec* const end);
void timespec_duration(const struct timespec* const begin, double* eval);
 
#ifdef __cplusplus
}
#endif

#endif
