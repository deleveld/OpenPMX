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

#ifndef VECTOR_C
#define VECTOR_C

#include <assert.h>
#include <string.h>
#include "vector.h"

void _vector_free(VECTOR_DATA * vect)
{
	if (vect->_capacity >= 0)
		free(vect->_data);
	memset(vect, 0, sizeof(*vect));
}

void _vector_reserve(VECTOR_DATA * vect, const int newsize, const size_t sizeofdata)
{
	assert(vect);
	assert(sizeofdata > 0);

	if (newsize > abs(vect->_capacity)) {

		/* fail reallocate when we are suppose to be fixed size buffer */
		assert(vect->_capacity >= 0);

		/* fail if it didnt work */
		vect->_data = realloc(vect->_data, newsize * sizeofdata);
		vect->_capacity = newsize;
		assert(vect->_data);
	}
}

void _vector_resize(VECTOR_DATA * vect, const int newsize, const size_t sizeofdata)
{
	assert(vect);
	assert(sizeofdata);

	if (newsize > abs(vect->_capacity)) {
		const int n = (newsize * 3) / 2;
		const int cap = (n > 16) ? n : 16;

		/* try to allocate some extra room */
		/* fail if it didnt work */
		vect->_data = realloc(vect->_data, cap * sizeofdata);
		vect->_capacity = cap;
		assert(vect->_data);
	}
	vect->_size = newsize;
}

void _vector_remove(VECTOR_DATA * vect, const int dindex, int num, const size_t sizeofdata)
{
	assert(vect);
	assert(sizeofdata);
	assert(dindex < vect->_size);

	int bytes;
	unsigned char *begin = (unsigned char *) vect->_data + dindex * sizeofdata;

	if ((dindex + num) > vect->_size)
		num = vect->_size - dindex;

	bytes = (vect->_size - dindex - num) * sizeofdata;
	if (bytes)
		memmove(begin, begin + num * sizeofdata, bytes);

	vect->_size -= num;
}

#endif
