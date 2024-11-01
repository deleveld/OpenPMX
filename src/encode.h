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

#ifndef OPENPMX_ENCODE_H
#define OPENPMX_ENCODE_H

#include "popmodel.h"
#include "omegainfo.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	POPMODEL popmodel;
	OMEGAINFO omegainfo;
	const int nparam;
	double offset[OPENPMX_THETA_MAX];
} ENCODE;

ENCODE encode_init(const POPMODEL* const popmodel);
void encode_popmodel(const POPMODEL* popmodel, ENCODE* const encode);
void encode_update_popmodel(ENCODE* encode, const double* x);

#ifdef __cplusplus
}
#endif

#endif
