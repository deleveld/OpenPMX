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

#ifndef OPENPMX_GETDELIM_H
#define OPENPMX_GETDELIM_H

#include <stdio.h>

#include "vector.h"

#ifdef __cplusplus
extern "C" {
#endif

ssize_t openpmx_getdelim(char **lineptr, size_t *n, int delim, FILE *stream);

void strip_firstlast_space(char* s);
 
typedef VECTOR(char*) VECPTR;

typedef enum {
	GET_DELIM_SEP_COMMA,		/* force comma separator */
	GET_DELIM_SEP_WHITESPACE,	/* force whitespace separator */
	GET_DELIM_SEP_ANY,			/* use comma separator if comma exists */
} GET_DELIM_SEP;

int get_delim_tokens(char* line, VECPTR* namevec, const GET_DELIM_SEP sep);

#ifdef __cplusplus
}
#endif

#endif
