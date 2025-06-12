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

#ifndef C22_H
#define C22_H

/* 22nd Century C https://gist.github.com/gheoan/9968496959e5c9cf27cb98b61679bd63 */
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

/* automatic type */
#ifdef __cplusplus
#define __auto_type auto
#endif
#define var __auto_type
#define let const __auto_type

/* string comparison */
#define streq(x, y) (strcmp(x, y) == 0)

/* safe countof() */
#ifdef __cplusplus
/* http://www.reedbeta.com/blog/cpp-compile-time-array-size/ with small modifications */
#define countof(x) (sizeof(countof_helper(x)))
template <typename T, size_t N>
char(&countof_helper(T(&)[N]))[N];
#else
/* https://stackoverflow.com/questions/19452971/array-size-macro-that-rejects-pointers */
#define countof(A) \
	_Generic(&(A), \
		typeof((A)[0]) **: (void)0, \
		default: sizeof(A) / sizeof((A)[0]))
#endif

/* looping */
#define forcount(index, count) for (int index=0, size ## __LINE__ = count; index < size ## __LINE__; ++index)
#define forarray(index, array) forcount(index, countof(array))

/* typesafe malloc and calloc */
#define mallocvar(t, n)	((t*)malloc((size_t)(n) * sizeof(t)))
#define callocvar(t, n) ((t*)calloc((size_t)(n), sizeof(t)))

/* min and max */
/* https://stackoverflow.com/questions/5595593/what-is-the-function-of-void-min1-min2-in-the-min-macro-in-kernel-h */
#ifndef min
#define min(x, y) ({                \
    typeof(x) _min1 = (x);          \
    typeof(y) _min2 = (y);          \
    (void) (&_min1 == &_min2);      \
    _min1 < _min2 ? _min1 : _min2; })
#endif

#ifndef max
#define max(x, y) ({                \
    typeof(x) _min1 = (x);          \
    typeof(y) _min2 = (y);          \
    (void) (&_min1 == &_min2);      \
    _min1 > _min2 ? _min1 : _min2; })
#endif

#endif
