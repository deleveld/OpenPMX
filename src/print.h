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

#ifndef OPENPMX_PRINT_H
#define OPENPMX_PRINT_H

#include <stdio.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

void print_serialize(const bool serial);

void openpmx_printf(FILE* stream1, FILE* stream2, const char* prefix, const char* format, ... );

void fatal(FILE* stream, const char* format, ... );
void warning(FILE* stream, const char* format, ... );
void info(FILE* stream, const char* format, ... );
 
#ifdef __cplusplus
}
#endif

#endif
