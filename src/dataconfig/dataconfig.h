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

#ifndef OPENPMX_DATACONFIG_INTERNAL_H
#define OPENPMX_DATACONFIG_INTERNAL_H

#include <math.h>

#include "openpmx.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct DATAINFO {
	const RECORD* record;
	double ID;
	double TIME;
	double DV;
	double EVID;
	double AMT;
	double RATE;
	double CMT;
	int CMT0;
	double DVLOW;
} DATAINFO;

typedef struct {
	const DATACONFIG* dataconfig;
	const int offsetID;
	const int offsetTIME;
	const int offsetDV;
	const int offsetEVID;
	const int offsetMDV;
	const int offsetAMT;
	const int offsetRATE;
	const int offsetCMT;
	const int offsetDVLOW;
	const bool _offset1;
	const int ndata;
	const int nindivid;
	const int nobs;
	const DATAINFO* datainfo;
} RECORDINFO;

RECORDINFO recordinfo_alloc(const DATACONFIG* const dataconfig);
void recordinfo_free(const RECORDINFO* recordinfo);

int structinfo_find_offset(const char* name, const STRUCTINFO* const structinfo);

static inline double DATA_FIELD(const void* p, const int offset)
{
	return *(const double*)((const char*)p + offset);
}

#ifdef __cplusplus
}
#endif

#endif
