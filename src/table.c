
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <ctype.h>
#include <dirent.h>

#include "openpmx.h"
#include "openpmx_internal.h"
#include "print.h"
#include "defines.h"
#include "advan/advan.h"
#include "utils/c22.h"
#include "utils/vector.h"
#include "utils/various.h"

#include "openpmx_compile_options.h"

static void strip_firstlast_space(char* s)
{
	let len = strlen(s);
	if (len == 0)
		return;

	/* remove space at end */
	char* a = s + len - 1;
	while (a != s && *a && isspace((int)*a)) {
		*a = 0;
		--a;
	}

	/* remove space at begin */
	a = s;
	while (*a && isspace((int)*a))
		++a;

	/* copy data over with classic strcpy
	 * which allows overlap */
	while ((*s++ = *a++)) { }
}

typedef struct {
	const IDATA* idata;
	const ADVANFUNCS* advanfuncs;
	const POPMODEL popmodel;
	const TABLECONFIG* tableconfig;
	char filename[PATH_MAX];
	char ext[NAME_MAX];

	const double zeroerr[OPENPMX_SIGMA_MAX];

	int printrow;
	int individ;
	int individ_i;
	FILE* stream;

	char* fields;						/* copy of the fields string */
	VECTOR(const char*) fieldnames;		/* pointer to each field */

	/* values for user */
	/* actually we dont have to calculate all of them at each table row
	 * we could do that when asked actually */
	const IMODEL* imodel;
	const PREDICTVARS* predictvars;
	const RECORD* record;
	const double* state;
	const double* eta;
	double yhat;
	double yhatvar;
	double pred;
	double obj;
	const double* err;
} TABLE;

static TABLE table_open(const IDATA* const idata,
						const ADVANFUNCS* const advanfuncs,
						const POPMODEL* const popmodel,
						const char* pmx_filename,
						int* tablenum,
						const char* _fields,
						const TABLECONFIG* const tableconfig)
{
	const char* filename = 0;
	char ext[PATH_MAX] = "";

	/* determine the file to open */
	if (tableconfig) {

		/* write to user define filename */
		if (tableconfig->filename) {
			if (tableconfig->name)
				fatal(0, "Table has both filename and name\n");
			filename = tableconfig->filename;

		/* with name, you prepend the control file name */
		} else if (tableconfig->name) {
			filename = pmx_filename;
			sprintf(ext, ".%s" OPENPMX_TABLEFILE, tableconfig->name);
		}
	}

	/* if no config then number the tables */
	if (!filename) {
		*tablenum += 1;
		filename = pmx_filename;
		sprintf(ext, ".table.%i" OPENPMX_TABLEFILE, *tablenum);
	}

	FILE* stream = stdout;
	if (filename)
		stream = results_fopen(filename, ext, "w");
	assert(stream);

	var ret = (TABLE) {
		.idata = idata,
		.advanfuncs = advanfuncs,
		.popmodel = *popmodel,
		.tableconfig = tableconfig,
		.filename = "",

		.zeroerr = { 0 },

		.printrow = 0,
		.individ = 0,
		.individ_i = 0,
		.stream = stream,
		.fields = strdup(_fields), /* make a copy for ourselves, when we tokenize our fields will be char* within this memory */
		/* fieldnames initialized after */

		.imodel = 0,
		.predictvars = 0,
		.record = 0,
		.state = 0,
		.eta = 0,
		.pred = 0,
		.yhat = 0,
		.err = 0,
	};
	/* save the full filename and ext */
	strcpy(ret.filename, filename);
	strcpy(ret.ext, ext);
	
	/* parse fields string and write the header */
	/* we have to tokenize a copy of the fields */
	char* token;
	char* rest = ret.fields;
	while ((token = strtok_r(rest, ", ", &rest))) {
		strip_firstlast_space(token);
		vector_append(ret.fieldnames, token);
		fprintf(stream, OPENPMX_HEADER_FORMAT, token);
	}
	fputc('\n', stream);

	return ret;
}

static void table_close(TABLE* const table)
{
	/* close the stream */
	if (table->stream)
		fclose(table->stream);
	table->stream = 0;

	/* cleanup */
	vector_free(table->fieldnames);
	free(table->fields);
}

/* table generation is single threaded because no integration takes place,
   only the saved state is used at each step so it should be I/O limited
   and not cpu limited. It would be very hard to make it multithreaded
   because the rows have to be output in the correct order */
static int table_row(TABLE* const table)
{
	assert(table);
	let idata = table->idata;
	let advanfuncs = table->advanfuncs;

	/* for first iteration setup, everything is 0 */
	if (table->record == 0) {
		table->individ = 0;
		table->individ_i = 0;
		table->printrow = 0;

	/* we are already handling a normal row */
	} else {
		bool firstonly = false;
		if (table->tableconfig)
			firstonly = table->tableconfig->firstonly;
		let nrecord = idata->individ[table->individ].nrecord;

		/* if firstonly then skip to the next individual */
		if (firstonly) {
			table->individ += 1;
			table->individ_i = 0;

		/* not firstonly but we are at the end of the indiviual */
		} else if (table->individ_i >= nrecord - 1) {
			/* not the last individual, so go to the next one */
			table->individ += 1;
			table->individ_i = 0;

		/* within an indivual we just go to the next row */
		} else {
			table->individ_i += 1;
		}
		/* stop at the end of the last individual, reset so we start over next time */
		if (table->individ >= idata->nindivid) {
			table->record = 0;
			table->individ = 0;
			table->individ_i = 0;
			table->printrow = 0;
			return 0;
		}

		/* go to the next row */
		table->printrow +=1;
	}

	/* prepare information about the table row for user */
	let record_size = advanfuncs->recordinfo.dataconfig->recordfields.size;
	let individ = &idata->individ[table->individ];
	table->imodel = (const IMODEL*)((char*)individ->imodel + table->individ_i * idata->imodel_size);
	table->predictvars = (const PREDICTVARS*)((char*)individ->predictvars + table->individ_i * idata->predictvars_size);
	table->record = RECORD_INDEX(individ->record, record_size, table->individ_i);
	table->state = individ->istate + table->individ_i * idata->nstate;
	table->eta = individ->eta;
	table->yhat = individ->yhat[table->individ_i];
	table->yhatvar = individ->yhatvar[table->individ_i];
	table->pred = individ->pred[table->individ_i];
	table->obj = individ->iobjfn;

	table->err = table->zeroerr;
	if (individ->isimerr) {
		let v = &individ->isimerr[table->individ_i * idata->nsigma];
		table->err = v;
	}

	return 1;
}

static double table_value(const TABLE* const table, const char* const name, const bool _offset1)
{
	let idata = table->idata;
	let individ = &idata->individ[table->individ];
	let advanfuncs = table->advanfuncs;

	let recordinfo = &advanfuncs->recordinfo;
	let recordfields = &recordinfo->dataconfig->recordfields;
	var recordoffset = structinfo_find_offset(name, recordfields);
	if (recordoffset >= 0)
		return DATA_FIELD(table->record, recordoffset);

	let modelfields = &advanfuncs->advanconfig->imodelfields;
	let modeloffset = structinfo_find_offset(name, modelfields);
	if (modeloffset >= 0)
		return DATA_FIELD(table->imodel, modeloffset);

	let predictfields = &advanfuncs->advanconfig->predictfields;
	let predictvarsoffset = structinfo_find_offset(name, predictfields);
	if (predictvarsoffset >= 0)
		return DATA_FIELD(table->predictvars, predictvarsoffset);

	if (strcmp(name, "yhat") == 0 ||
		strcmp(name, "YHAT") == 0 ||
		strcmp(name, "y") == 0 ||
		strcmp(name, "Y") == 0)
		return table->yhat;
	if (strcmp(name, "yhatvar") == 0 ||
		strcmp(name, "YHATVAR") == 0)
		return table->yhatvar;
	if (strcmp(name, "PRED") == 0 ||
		strcmp(name, "pred") == 0)
		return table->pred;
	if (strcmp(name, "obj") == 0 ||
		strcmp(name, "OBJ") == 0)
		return table->obj;
	if (strcmp(name, "ineval") == 0 ||
		strcmp(name, "INEVAL") == 0)
		return individ->ineval;
	if (strcmp(name, "cwres") == 0 ||
		strcmp(name, "CWRES") == 0) {
		let dv = RECORDINFO_DV(recordinfo, table->record);
		return (table->yhat - dv) / sqrt(table->yhatvar);
	}

	if (strcmp(name, "evid") == 0 ||
		strcmp(name, "EVID") == 0) 
		return RECORDINFO_EVID(recordinfo, table->record);
	if (strcmp(name, "MDV") == 0)
		return RECORDINFO_MDV(recordinfo, table->record);
	if (strcmp(name, "EVID") == 0)
		return RECORDINFO_EVID(recordinfo, table->record);
	if (strcmp(name, "AMT") == 0)
		return RECORDINFO_AMT(recordinfo, table->record);
	if (strcmp(name, "RATE") == 0)
		return RECORDINFO_RATE(recordinfo, table->record);
	if (strcmp(name, "CMT") == 0)
		return RECORDINFO_CMT(recordinfo, table->record);
		
	var i = -1;

	/* normal C access */
	if (_offset1 == false) {
		if (sscanf(name, "theta[%i]", &i)) {
			assert(i >= 0 && i < idata->ntheta);
			return table->popmodel.theta[i];
		}
		if (sscanf(name, "eta[%i]", &i)) {
			assert(i >= 0 && i < idata->nomega);
			return table->eta[i];
		}
		if ((sscanf(name, "a(%i)", &i) == 1) ||
			(sscanf(name, "state(%i)", &i) == 1)) {
			assert(i >= 0 && i < idata->nstate);
			return table->state[i];
		}
		if (sscanf(name, "err[%i]", &i) == 1) {
			assert(i >= 0 && i < idata->nsigma);
			return table->err[i];
		}

	/* NONMEM like access with 1-offset vectors */
	} else {
		if (sscanf(name, "THETA(%i)", &i) == 1) {
			assert(i-1 >= 0 && i-1 < idata->ntheta);
			return table->popmodel.theta[i-1];
		}
		if (sscanf(name, "ETA(%i)", &i) == 1) {
			assert(i-1 >= 0 && i-1 < idata->nomega);
			return table->eta[i-1];
		}
		if ((sscanf(name, "A(%i)", &i) == 1) ||
			(sscanf(name, "STATE(%i)", &i) == 1)) {
			assert(i-1 >= 0 && i-1 < advanfuncs->nstate);
			return table->state[i-1];
		}
		if (sscanf(name, "ERR(%i)", &i) == 1) {
			assert(i-1 >= 0 && i-1 < idata->nsigma);
			return table->err[i-1];
		}
	}
	
	fatal(0, "Cannot find table field \"%s\"\n", name);
	return NAN;
}

void pmx_table(OPENPMX* pmx,
			   const char* fields,
			   const TABLECONFIG* const tableconfig)
{
	pmxstate_ensure(pmx);
	
	let popmodel = popmodel_init(pmx);
	TABLE table = table_open(&pmx->state->idata,
							 pmx->state->advanfuncs,
							 &popmodel,
							 pmx->filename,
							 &pmx->state->tablecount,
							 fields,
							 tableconfig);

	FILE* f = table.stream;
	if (f == 0)
		f = stdout;

	/* write out fields */
	while (table_row(&table)) {
		forvector(i, table.fieldnames) {
			let name = table.fieldnames.data[i];
			let val = table_value(&table, name, pmx->_offset1);
			fprintf(f, OPENPMX_TABLE_FORMAT, val);
		}
		fputc('\n', f);
	}
	table_close(&table);
}

