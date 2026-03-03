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
 
/// This file initializes a RECORDINFO object from a DATACONFIG object.
/// In this way the offsets of some important RECORD fields are cached
/// for rapid access, 

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

#include "dataconfig/dataconfig.h"
#include "utils/c22.h"
#include "utils/vector.h"

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

static inline double RECORDINFO_AMT(const RECORDINFO* const recordinfo, const RECORD* p)
{
	if (recordinfo->offsetAMT != -1)
		return *(const double*)((const char*)p + recordinfo->offsetAMT);
	return 0.;
}

static inline double RECORDINFO_RATE(const RECORDINFO* const recordinfo, const RECORD* p)
{
	if (recordinfo->offsetRATE != -1)
		return *(const double*)((const char*)p + recordinfo->offsetRATE);
	return 0.;
}

static inline int RECORDINFO_MDV(const RECORDINFO* const recordinfo, const RECORD* p)
{
	if (recordinfo->offsetMDV != -1)
		return (int)(*(const double*)((const char*)p + recordinfo->offsetMDV));
	if (isnan(RECORDINFO_DV(recordinfo, p)))
		return 1;
	return 0;
}

static inline int RECORDINFO_EVID(const RECORDINFO* const recordinfo, const RECORD* p)
{
	if (recordinfo->offsetEVID != -1) 
		return (int)(*(const double*)((const char*)p + recordinfo->offsetEVID));
	if (RECORDINFO_MDV(recordinfo, p) == 0)
		return 0;
	if (RECORDINFO_AMT(recordinfo, p) != 0.)
		return 1;
	return 2;
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
		return 0;

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

RECORDINFO recordinfo_alloc(const DATACONFIG* const dataconfig)
{
    let data = dataconfig->records;
    let recordfields = &dataconfig->recordfields;

    /* Bootstrap an incomplete struct to utilize accessor functions for counting */
    let temp = (RECORDINFO) {
        .dataconfig     = dataconfig,
        .offsetID       = structinfo_find_offset("ID", recordfields),
        .offsetTIME     = structinfo_find_offset("TIME", recordfields),
        .offsetDV       = structinfo_find_offset("DV", recordfields),
        .offsetEVID     = structinfo_find_offset("EVID", recordfields),
        .offsetMDV      = structinfo_find_offset("MDV", recordfields),
        .offsetAMT      = structinfo_find_offset("AMT", recordfields),
        .offsetRATE     = structinfo_find_offset("RATE", recordfields),
        .offsetCMT      = structinfo_find_offset("CMT", recordfields),
        .offsetDVLOW    = structinfo_find_offset("DVLOW", recordfields),
        ._offset1       = dataconfig->_offset1,
        .ndata          = 0,
        .nindivid       = 0,
        .nobs           = 0,
    };

    let nrecords = dataconfig->nrecords;
    var i = 0;
    var ndata = 0;
    var nobs = 0;
    var nindivid = 0;
    VECTOR(DATAINFO) datainfo = { 0 };

    while (i < nrecords) {
        let first_rec_in_indiv = RECORDINFO_INDEX(&temp, data, i);
        let thisid = RECORDINFO_ID(&temp, first_rec_in_indiv);

        /* check for NaN ID at the start of an individual block */
        if (isnan(thisid)) {
            fprintf(stderr, "error: NaN ID detected at record index %d.\n", i);
            exit(EXIT_FAILURE);
        }

        int nidata = 0;
        while (i + nidata < nrecords) {
            let curr_rec = RECORDINFO_INDEX(&temp, data, i + nidata);
            let current_id = RECORDINFO_ID(&temp, curr_rec);

            /* check for NaN ID within a data block */
            if (isnan(current_id)) {
                fprintf(stderr, "error: NaN ID detected within data at record index %d.\n", i + nidata);
                exit(EXIT_FAILURE);
            }

            /* boundary check for new individual */
            if (current_id != thisid) 
                break;

            /* Count observations: Non-NaN DV and EVID 0 */
            let dv = RECORDINFO_DV(&temp, curr_rec);
            let evid = RECORDINFO_EVID(&temp, curr_rec);
            if (!isnan(dv) && evid == 0) 
                ++nobs;

            vector_append(datainfo,
				(DATAINFO) {
					.record = curr_rec,
					.ID 	= RECORDINFO_ID(&temp, curr_rec),
					.TIME 	= RECORDINFO_TIME(&temp, curr_rec),
					.DV 	= RECORDINFO_DV(&temp, curr_rec),
					.EVID 	= RECORDINFO_EVID(&temp, curr_rec),
					.AMT 	= RECORDINFO_AMT(&temp, curr_rec),
					.RATE 	= RECORDINFO_RATE(&temp, curr_rec),
					.CMT 	= RECORDINFO_CMT(&temp, curr_rec),
					.CMT0 	= RECORDINFO_CMT_0offset(&temp, curr_rec),
					.DVLOW 	= RECORDINFO_DVLOW(&temp, curr_rec),
				}
			);

            ++nidata;
            ++ndata;
        }

        ++nindivid;
        i += nidata;
    }

    return (RECORDINFO) {
        .dataconfig  = temp.dataconfig,
        .offsetID    = temp.offsetID,
        .offsetTIME  = temp.offsetTIME,
        .offsetDV    = temp.offsetDV,
        .offsetEVID  = temp.offsetEVID,
        .offsetMDV   = temp.offsetMDV,
        .offsetAMT   = temp.offsetAMT,
        .offsetRATE  = temp.offsetRATE,
        .offsetCMT   = temp.offsetCMT,
        .offsetDVLOW = temp.offsetDVLOW,
        ._offset1    = temp._offset1,
        .ndata       = ndata,
        .nindivid    = nindivid,
        .nobs        = nobs,
        .datainfo    = datainfo.ptr,
    };
}

void recordinfo_free(const RECORDINFO* recordinfo)
{
	free((void*)recordinfo->datainfo);
}

int structinfo_find_offset(const char* name, const STRUCTINFO* const structinfo)
{
	for (int i=0; i<OPENPMX_FIELDS_MAX; i++) {
		if (structinfo->field[i].name[0] == '\0')
			return -1;
		if (strcasecmp(name, structinfo->field[i].name) == 0)
			return structinfo->field[i].offset;
	}
	return -1;
}
