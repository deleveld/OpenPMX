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

#define VECTOR(type) \
	union { \
		VECTOR_DATA _vector; \
		type* data; \
	}

#define vector_free(mxvect) _vector_free(&(mxvect)._vector);

#define vector_buffer(mxvect, _a, _b)				\
	do {											\
		VECTOR_DATA* _mxvect = &(mxvect)._vector; 	\
		assert(_mxvect->_data == 0);				\
		assert(_mxvect->_size == 0);				\
		assert(_mxvect->_capacity == 0);			\
		_mxvect->_data = &(_a)[0];					\
		_mxvect->_capacity = -(_b);					\
		assert(_mxvect->_capacity < 0);				\
	} while (0)

#define vector_size(mxvect)			((const int)((mxvect)._vector._size))
#define vector_capacity(mxvect)		((const int)((mxvect)._vector._capacity))

#define vector_at(mxvect, ...)		((mxvect).data[__VA_ARGS__])

#define vector_first(mxvect)	 	((mxvect).data[0])
#define vector_last(mxvect) 		((mxvect).data[vector_size(mxvect)-1])

#define vector_resize(mxvect, n)	_vector_resize(&(mxvect)._vector, (n), sizeof((mxvect).data[0]))
#define vector_reserve(mxvect, n)	_vector_reserve(&(mxvect)._vector, (n), sizeof((mxvect).data[0]))

#define vector_append(mxvect, ...) 							\
	do { 													\
		const int vector_append_n1 = vector_size(mxvect); 	\
		vector_resize(mxvect, vector_append_n1 + 1); 		\
		(mxvect).data[vector_append_n1] = (__VA_ARGS__); 	\
	} while (0)

#define vector_appendn(mxvect, dat, n) 						\
	do { 													\
		const int mx__n = (int)(n); 						\
		const int mx__n1 = vector_size(mxvect) + mx__n; 	\
		vector_resize(mxvect, mx__n1); 						\
		memcpy(&(mxvect).data[mx__n1 - mx__n], (dat), mx__n * sizeof((mxvect).data[0])); \
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
