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
	const int offsetDVLOW;
	const bool _offset1;
	const int ndata;
	const int nindivid;
	const int nobs;
} RECORDINFO;

RECORDINFO recordinfo_init(const DATACONFIG* const dataconfig);

static inline double DATA_FIELD(const void* p, const int offset)
{
	return *(const double*)((const char*)p + offset);
}

static inline const RECORD* RECORD_INDEX(const RECORD* p, const int size, const int i)
{
	return (const RECORD*)((const char*)p + size * (int)i);
}

static inline double RECORDINFO_ID(const RECORDINFO* const recordinfo, const RECORD* p)
{
	return *(const double*)((const char*)p + recordinfo->offsetID);
}

static inline double RECORDINFO_TIME(const RECORDINFO* const recordinfo, const RECORD* p)
{
	return *(const double*)((const char*)p + recordinfo->offsetTIME);
}

static inline double RECORDINFO_DV(const RECORDINFO* const recordinfo, const RECORD* p)
{
	return *(const double*)((const char*)p + recordinfo->offsetDV);
}

static inline int RECORDINFO_MDV(const RECORDINFO* const recordinfo, const RECORD* p)
{
	if (recordinfo->offsetMDV == -1) {
		return (isnan(RECORDINFO_DV(recordinfo, p)) ? 1 : 0);
	} else
		return (int)(*(const double*)((const char*)p + recordinfo->offsetMDV));
}

static inline int RECORDINFO_EVID(const RECORDINFO* const recordinfo, const RECORD* p)
{
	if (recordinfo->offsetEVID == -1)
		return (RECORDINFO_MDV(recordinfo, p) != 0 ? 1 : 0);
	else
		return (int)(*(const double*)((const char*)p + recordinfo->offsetEVID));
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

static inline int RECORDINFO_CMT_0offset(const RECORDINFO* const recordinfo, const RECORD* p)
{
	if (recordinfo->offsetCMT == -1)
		return 0.;

	const int ret = (int)(*(const double*)((const char*)p + recordinfo->offsetCMT));
	return (recordinfo->_offset1) ? (ret - 1) : ret;
}

static inline double RECORDINFO_CMT(const RECORDINFO* const recordinfo, const RECORD* p)
{
	if (recordinfo->offsetCMT == -1)
		return 0.;

	return *(const double*)((const char*)p + recordinfo->offsetCMT);
}

static inline double RECORDINFO_DVLOW(const RECORDINFO* const recordinfo, const RECORD* p)
{
	if (recordinfo->offsetDVLOW == -1)
		return 0.;

	return *(const double*)((const char*)p + recordinfo->offsetDVLOW);
}

static inline const RECORD* RECORDINFO_INDEX(const RECORDINFO* const recordinfo, const RECORD* p, const int i)
{
	return RECORD_INDEX(p, recordinfo->dataconfig->recordfields.size, i);
}

#ifdef __cplusplus
}
#endif

#endif
