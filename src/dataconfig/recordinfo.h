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

#ifndef OPENPMX_RECORDINFO_H
#define OPENPMX_RECORDINFO_H

#include "dataconfig/dataconfig.h"

#include <math.h>
#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

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
	const bool _offset1;
	const int ndata;
	const int nindivid;
	const int nobs;
} RECORDINFO;

RECORDINFO recordinfo_init(const DATACONFIG* const dataconfig);

inline double DATA_FIELD(const void* p, const int offset)
{
	return *(const double*)((const char*)p + offset);
}

inline const RECORD* RECORD_INDEX(const RECORD* p, const int size, const int i)
{
	return (const RECORD*)((const char*)p + size * (int)i);
}

inline double RECORDINFO_ID(const RECORDINFO* const recordinfo, const RECORD* p)
{
	return *(const double*)((const char*)p + recordinfo->offsetID);
}

inline double RECORDINFO_TIME(const RECORDINFO* const recordinfo, const RECORD* p)
{
	return *(const double*)((const char*)p + recordinfo->offsetTIME);
}

inline double RECORDINFO_DV(const RECORDINFO* const recordinfo, const RECORD* p)
{
	return *(const double*)((const char*)p + recordinfo->offsetDV);
}

inline double RECORDINFO_MDV(const RECORDINFO* const recordinfo, const RECORD* p)
{
	if (recordinfo->offsetMDV == -1) {
		return (isnan(RECORDINFO_DV(recordinfo, p)) ? 1. : 0.);
	} else
		return *(const double*)((const char*)p + recordinfo->offsetMDV);
}

inline double RECORDINFO_EVID(const RECORDINFO* const recordinfo, const RECORD* p)
{
	if (recordinfo->offsetEVID == -1)
		return (RECORDINFO_MDV(recordinfo, p) != 0 ? 1. : 0.);
	else
		return *(const double*)((const char*)p + recordinfo->offsetEVID);
}

static inline double RECORDINFO_AMT(const RECORDINFO* const recordinfo, const RECORD* p)
{
	if (recordinfo->offsetAMT == -1)
		return 0.;
	else
		return *(const double*)((const char*)p + recordinfo->offsetAMT);
}

static inline double RECORDINFO_RATE(const RECORDINFO* const recordinfo, const RECORD* p)
{
	if (recordinfo->offsetRATE == -1)
		return 0.;
	else
		return *(const double*)((const char*)p + recordinfo->offsetRATE);
}

static inline double RECORDINFO_CMT(const RECORDINFO* const recordinfo, const RECORD* p)
{
	if (recordinfo->offsetCMT == -1)
		return 0.;

	const double ret = *(const double*)((const char*)p + recordinfo->offsetCMT);
	if (recordinfo->_offset1)
		return ret - 1;
	return ret;
}

inline const RECORD* RECORDINFO_INDEX(const RECORDINFO* const recordinfo, const RECORD* p, const int i)
{
	return RECORD_INDEX(p, recordinfo->dataconfig->recordfields.size, i);
}

#ifdef __cplusplus
}
#endif

#endif
