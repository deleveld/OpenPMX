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
		type* mutptr;				\
	}

#define vector_alloc(_mxvect) 		{  }
#define vector_free(_mxvect) 		_vector_free(&(_mxvect)._vector);

#define vector_resize(_mxvect, _n)	_vector_resize(&(_mxvect)._vector, (_n), sizeof((_mxvect).ptr[0]))
#define vector_reserve(_mxvect, _n)	_vector_reserve(&(_mxvect)._vector, (_n), sizeof((_mxvect).ptr[0]))

#define vector_append(_mxvect, ...)                                    \
	do {                                                               \
		__auto_type _mxvectptr = &(_mxvect);                    \
		__auto_type mx__n = _mxvectptr->size;                            \
		vector_resize(*_mxvectptr, mx__n + 1);                         \
		_mxvectptr->mutptr[mx__n] = (__VA_ARGS__);                     \
	} while (0)

#define vector_appendn(_mxvect, _dat, _n)                              \
	do {                                                               \
		__auto_type _mxvectptr = &(_mxvect);                    \
		__auto_type _datptr = (_dat);                                 \
		__auto_type _datn = (int)(_n);                                   \
		__auto_type _oldn = _mxvectptr->size;                            \
		vector_resize(*_mxvectptr, _oldn + _datn);                     \
		memcpy(&_mxvectptr->mutptr[_oldn], _datptr, _datn * sizeof(_mxvectptr->ptr[0])); \
	} while (0)

#define vector_remove(_mxvect, _ind, _n)                               \
	_vector_remove(&(_mxvect)._vector, (_ind), (_n), sizeof((_mxvect).ptr[0]))

#define vector_remove_first(_mxvect) vector_remove((_mxvect), 0, 1)
#define vector_remove_last(_mxvect) vector_resize((_mxvect), (_mxvect).size - 1)

#define vector_insert(_mxvect, _ind, ...)                              \
    do {                                                               \
        __auto_type _mxvectptr = &(_mxvect);                           \
        __auto_type _mxind = (_ind);                                   \
        _vector_insert(&_mxvectptr->_vector, _mxind, 1,                \
                       sizeof(_mxvectptr->ptr[0]));                    \
        _mxvectptr->mutptr[_mxind] = (__VA_ARGS__);                    \
    } while (0)
    
#define vector_insertn(_mxvect, _ind, _dat, _n)                        \
    do {                                                               \
        __auto_type _mxvectptr = &(_mxvect);                           \
        __auto_type _datptr = (_dat);                                  \
        __auto_type _mxind = (_ind);                                   \
		__auto_type _datn = (int)(_n);                                 \
        _vector_insert(&_mxvectptr->_vector, _mxind, _datn,            \
                       sizeof(_mxvectptr->ptr[0]));                    \
        memcpy(&_mxvectptr->mutptr[_mxind], _datptr,                   \
               _datn * sizeof(_mxvectptr->ptr[0]));                    \
    } while (0)
    
#define forvector(_index, _mxvect)                          \
    for (__auto_type _vptr = &(_mxvect);                    \
         _vptr;                                             \
         _vptr = NULL)                                      \
        for (int _index = 0; _index < _vptr->size; _index++)

#define forvector_val(_val, _mxvect)                        \
    for (__auto_type _vptr = &(_mxvect);                    \
         _vptr;                                             \
         _vptr = NULL)                                      \
        for (__auto_type _it = _vptr->ptr;                  \
             _it != _vptr->ptr + _vptr->size;               \
             ++_it)                                         \
            for (int _once = 1; _once; _once = 0)           \
                for (__auto_type _val = *_it; _once; _once = 0)

#define forvector_ptr(_ptr, _mxvect)                        \
    for (__auto_type _vptr = &(_mxvect);                    \
         _vptr;                                             \
         _vptr = NULL)                                      \
        for (__auto_type _ptr = _vptr->ptr;                 \
             _ptr != _vptr->ptr + _vptr->size;              \
             ++_ptr)

#define forvector_mutptr(_ptr, _mxvect)                     \
    for (__auto_type _vptr = &(_mxvect);                    \
         _vptr;                                             \
         _vptr = NULL)                                      \
        for (__auto_type _ptr = _vptr->mutptr;              \
             _ptr != _vptr->mutptr + _vptr->size;           \
             ++_ptr)      
             
	void _vector_free(VECTOR_DATA * vect);
	void _vector_resize(VECTOR_DATA * vect, const int newsize, const size_t sizeofdata);
	void _vector_reserve(VECTOR_DATA * vect, const int newsize, const size_t sizeofdata);
	void _vector_remove(VECTOR_DATA * vect, const int dindex, int num, const size_t sizeofdata);
	void _vector_insert(VECTOR_DATA * vect, const int dindex, const int num, const size_t sizeofdata);

#ifdef __cplusplus
}
#endif
#endif
