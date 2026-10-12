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

#ifndef OPENPMX_COVARIANCE_H
#define OPENPMX_COVARIANCE_H

#include "openpmx.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef typeof(((COVARIANCECONFIG*)0)->type) COVARIANCETYPE;

/* Result of pmx_covariance() stored in the OPENPMX state. 'limit' is an
 * owned array of nlimit confidence limits according to the change in
 * objective function. */
typedef struct COVARIANCE {
	COVARIANCETYPE type;
	const COVLIMIT* limit;
	int nlimit;
	double dobjfn;
} COVARIANCE;

/* Releases the limit array and resets it, safe to call repeatedly. */
void pmx_covariance_free(COVARIANCE* cov);

#ifdef __cplusplus
}
#endif

#endif
