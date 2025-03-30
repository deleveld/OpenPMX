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

#ifndef VECTOR_H
#define VECTOR_H

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

	typedef struct {
		void *_data;
		int _size;
		int _capacity;
	} VECTOR_DATA;

#define VECTOR(type)                                                   \
	union {                                                            \
		VECTOR_DATA _vector;                                           \
		type* data;                                                    \
	}

#define vector_alloc(mxvect) 		{  }
#define vector_init_buffer(mxvect, _a, _b)                             \
do {                                                                   \
	(mxvect)._vector = (VECTOR_DATA) { ._data=&((_a)[0]), ._size = 0, ._capacity=-abs(_b), };\
} while (0)

#define vector_free(mxvect) _vector_free(&(mxvect)._vector);

#define vector_size(mxvect)			((const int)((mxvect)._vector._size))
#define vector_capacity(mxvect)		((const int)abs((mxvect)._vector._capacity))

#define vector_at(mxvect, ...)		((mxvect).data[__VA_ARGS__])

#define vector_first(mxvect)	 	((mxvect).data[0])
#define vector_last(mxvect) 		((mxvect).data[vector_size(mxvect)-1])

#define vector_resize(mxvect, n)	_vector_resize(&(mxvect)._vector, (n), sizeof((mxvect).data[0]))
#define vector_reserve(mxvect, n)	_vector_reserve(&(mxvect)._vector, (n), sizeof((mxvect).data[0]))

#define vector_append(mxvect, ...)                                     \
	do {                                                               \
		const int mx__n = vector_size(mxvect);                         \
		vector_resize(mxvect, mx__n + 1);                              \
		(mxvect).data[mx__n] = (__VA_ARGS__);                          \
	} while (0)

#define vector_appendn(mxvect, _dat, _n)                               \
	do {                                                               \
		const int mx__n = vector_size(mxvect);                         \
		vector_resize(mxvect, mx__n + (int)_n);                        \
		for (int mx__i=0; mx__i<(int)_n; ++mx__i)                      \
			(mxvect).data[mx__n + mx__i] = (_dat)[mx__i];              \
	} while (0)

#define vector_remove(mxvect, ind, n) _vector_remove(&(mxvect)._vector, (ind), (n), (int)sizeof((mxvect).data[0]))

#define vector_remove_first(mxvect) vector_remove(mxvect, 0, 1)
#define vector_remove_last(mxvect) vector_resize(mxvect, vector_size(mxvect) - 1)

#define forvector(index, mxvect) for(int index=0; index < vector_size(mxvect); index++)

	void _vector_free(VECTOR_DATA * vect);
	void _vector_resize(VECTOR_DATA * vect, const int newsize, const size_t sizeofdata);
	void _vector_reserve(VECTOR_DATA * vect, const int newsize, const size_t sizeofdata);
	void _vector_remove(VECTOR_DATA * vect, const int dindex, int num, const size_t sizeofdata);

#ifdef __cplusplus
}
#endif
#endif
