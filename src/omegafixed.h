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

#ifndef OPENPMX_OMEGAFIXED_H
#define OPENPMX_OMEGAFIXED_H

#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned char OMEGAFIXED;
enum {
	OMEGAFIXED_ESTIMATE = 0,
	OMEGAFIXED_FIXED,
	OMEGAFIXED_SAME,
};

/* Conversion to the value in the .ext file from OMEGAFIXED. This is 
 * used when writing the .ext file header in popmodel.c */
static inline double omegafixed_to_ext_fixedval(const OMEGAFIXED f)
{
	int v = -1.; 
	if (f == OMEGAFIXED_ESTIMATE)
		v = 0.;
	else if (f == OMEGAFIXED_FIXED)
		v = 1.;
	else if (f == OMEGAFIXED_SAME)
		v = 2.;
	return v;
}

/* Conversion to OMEGAFIXED from the value in the .ext file. This is
 * used in reload.c */
static inline OMEGAFIXED omegafixed_from_ext_fixedval(const int val)
{
	OMEGAFIXED ret = OMEGAFIXED_ESTIMATE; /* for val == 0 */
	if (val == 1)
		ret = OMEGAFIXED_FIXED;
	else if (val == 2)
		ret = OMEGAFIXED_SAME;
	return ret;
}

#ifdef __cplusplus
}
#endif

#endif
