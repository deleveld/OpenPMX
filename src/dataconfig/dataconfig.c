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

RECORDINFO recordinfo_init(const DATACONFIG* const dataconfig)
{
	const RECORD* data = dataconfig->records;
	const STRUCTINFO* recordfields = &dataconfig->recordfields;

	/* construct an incomplete struct so we can use the assesor functions so that we can count the number
	 * of individuals and observations in the data.
	 * This is done in this indirect somewhat awkward way so that the RECORDINFO struct can have const members */
	const RECORDINFO temp = (RECORDINFO) {
		.dataconfig 	= dataconfig,
		.offsetID 		= structinfo_find_offset("ID", recordfields),
		.offsetTIME		= structinfo_find_offset("TIME", recordfields),
		.offsetDV		= structinfo_find_offset("DV", recordfields),
		.offsetEVID		= structinfo_find_offset("EVID", recordfields),
		.offsetMDV		= structinfo_find_offset("MDV", recordfields),
		.offsetAMT		= structinfo_find_offset("AMT", recordfields),
		.offsetRATE		= structinfo_find_offset("RATE", recordfields),
		.offsetCMT		= structinfo_find_offset("CMT", recordfields),
		.offsetDVLOW	= structinfo_find_offset("DVLOW", recordfields),
		._offset1		= dataconfig->_offset1,
		.ndata			= 0,
		.nindivid		= 0,
		.nobs			= 0,
	};
	
	/* count the number of individuals and observations using the incompletely constructed RECORDINFO */
	const int nrecords = dataconfig->nrecords;
	int i = 0;
	int ndata = 0;
	int nobs = 0;
	int nindivid = 0;
	while (i < nrecords) {
		int nidata = 0;
		double thisid = RECORDINFO_ID(&temp, RECORDINFO_INDEX(&temp, data, i));
		if (isnan(thisid))
			break;
		while (i+nidata < nrecords && RECORDINFO_ID(&temp, RECORDINFO_INDEX(&temp, data, i+nidata)) == thisid) {
			const RECORD* ri = RECORDINFO_INDEX(&temp, data, i+nidata);
			const double dv = RECORDINFO_DV(&temp, ri);
			const int evid = RECORDINFO_EVID(&temp, ri);
			if (!isnan(dv) && evid == 0)
				++nobs;
			++nidata;
			++ndata;
		}
		++nindivid;
		i += nidata;
	}

	return (RECORDINFO) {
		.dataconfig = temp.dataconfig,
		.offsetID 	= temp.offsetID,
		.offsetTIME = temp.offsetTIME,
		.offsetDV	= temp.offsetDV,
		.offsetEVID = temp.offsetEVID,
		.offsetMDV 	= temp.offsetMDV,
		.offsetAMT 	= temp.offsetAMT,
		.offsetRATE = temp.offsetRATE,
		.offsetCMT 	= temp.offsetCMT,
		.offsetDVLOW = temp.offsetDVLOW,
		._offset1 = temp._offset1,
		.ndata = ndata,
		.nindivid = nindivid,
		.nobs = nobs,
	};
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
