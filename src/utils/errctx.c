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

#include <string.h>
#include <stdio.h>

#include "errctx.h"

static void verrctx(ERRCTX *errctx, const char *format, va_list args)
{
	size_t capacity = sizeof(errctx->errmsg);
	if (errctx->len >= capacity - 1) 
		return;

	size_t available = capacity - errctx->len;
	int written = vsnprintf(&errctx->errmsg[errctx->len], available, format, args);
	if (written > 0) {
		if ((size_t)written >= available) 
			errctx->len = capacity - 1; 
		else 
			errctx->len += (size_t)written;
	}
}

void add_errctx(ERRCTX *errctx, const char *format, ...)
{
    va_list args;
    va_start(args, format);
    verrctx(errctx, format, args);
    va_end(args);
}
