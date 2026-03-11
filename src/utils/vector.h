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

#define VECTOR(type)				\
	union {							\
		VECTOR_DATA _vector;		\
		struct {					\
			const type* const ptr;	\
			const int size;			\
			const int capacity;		\
		};							\
		type* rawptr;				\
	}

#define vector_alloc(_mxvect) 		{  }
#define vector_free(_mxvect) 		_vector_free(&(_mxvect)._vector);

/* for a fixed buffer that does not realloc or need freeing */
#define vector_buffer(ptr, capacity)			{	._vector = _vector_buffer_init(0, capacity), .rawptr = ptr, }
#define vector_buffer_init(ptr, size, capacity)	{	._vector = _vector_buffer_init(size, capacity), .rawptr = ptr,  }
static inline VECTOR_DATA _vector_buffer_init(const int size, const int capacity)
{
	assert(size >= 0);
	assert(capacity > 0);
	assert(size <= capacity);
	return (VECTOR_DATA) {
		._size = size,
		._capacity = -capacity,
	};
}

#define vector_first(_mxvect)	 	((_mxvect).ptr[0])
#define vector_last(_mxvect) 		((_mxvect).ptr[(_mxvect).size-1])

#define vector_resize(_mxvect, _n)	_vector_resize(&(_mxvect)._vector, (_n), sizeof((_mxvect).ptr[0]))
#define vector_reserve(_mxvect, _n)	_vector_reserve(&(_mxvect)._vector, (_n), sizeof((_mxvect).ptr[0]))

#define vector_append(_mxvect, ...)                                    \
	do {                                                               \
		typeof(&(_mxvect)) _mxvectptr = &(_mxvect);                    \
		const int mx__n = _mxvectptr->size;                            \
		vector_resize(*_mxvectptr, mx__n + 1);                         \
		_mxvectptr->rawptr[mx__n] = (__VA_ARGS__);                     \
	} while (0)

#define vector_appendn(_mxvect, _dat, _n)                              \
	do {                                                               \
		typeof(&(_mxvect)) _mxvectptr = &(_mxvect);                    \
		typeof(_dat) _datptr = (_dat);                                 \
		const int _datn = (int)(_n);                                   \
		const int _oldn = _mxvectptr->size;                            \
		vector_resize(*_mxvectptr, _oldn + _datn);                     \
		for (int _i=0; _i<_datn; _i++)                                 \
			_mxvectptr->rawptr[_oldn + _i] = _datptr[_i];              \
	} while (0)

#define vector_remove(_mxvect, _ind, _n) _vector_remove(&(_mxvect)._vector, (_ind), (_n), (int)sizeof((_mxvect).ptr[0]))

#define vector_remove_first(_mxvect) vector_remove((_mxvect), 0, 1)
#define vector_remove_last(_mxvect) vector_resize((_mxvect), (_mxvect).size - 1)

#define forvector(_index, _mxvect) for(int _index=0; _index < (_mxvect).size; _index++)

	void _vector_free(VECTOR_DATA * vect);
	void _vector_resize(VECTOR_DATA * vect, const int newsize, const size_t sizeofdata);
	void _vector_reserve(VECTOR_DATA * vect, const int newsize, const size_t sizeofdata);
	void _vector_remove(VECTOR_DATA * vect, const int dindex, int num, const size_t sizeofdata);

#ifdef __cplusplus
}
#endif
#endif
