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

/// This file implements writing tables. 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <ctype.h>
#include <dirent.h>

#include "openpmx.h"
#include "pmxstate.h"
#include "print.h"
#include "defines.h"
#include "advan/advan.h"
#include "utils/c22.h"
#include "utils/vector.h"

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

	const double zeroerr[OPENPMX_SIGMA_MAX];

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

#ifndef NAME_MAX
#define NAME_MAX 256
#endif

static TABLE table_open(const IDATA* const idata,
						const ADVANFUNCS* const advanfuncs,
						const POPMODEL* const popmodel,
						const char* pmx_filename,
						int* tablecount,
						const char* _fields,
						const TABLECONFIG* const tableconfig)
{
	const char* filename = 0;
	char ext[NAME_MAX] = "";

	/* determine the file to open */
	if (tableconfig) {

/// If a filename is defined in the TABLECONFIG then the output is
/// written to that filename.  
		/* write to user define filename */
		if (tableconfig->filename) {
			if (tableconfig->name)
				fatal(0, "Table has both filename and name\n");
			if (tableconfig->stream)
				fatal(0, "Table has both filename and stream\n");
			filename = tableconfig->filename;

/// If a name is defined in the TABLECONFIG then the output is
/// written to a filename which is the name in the PMX object extended
/// with the given name.   
		} else if (tableconfig->name) {
			if (tableconfig->stream)
				fatal(0, "Table has both name and stream\n");
			filename = pmx_filename;
			if (!filename)
				filename = "";
			sprintf(ext, "%s", tableconfig->name);
		}
	}

	FILE* stream = 0;
/// if a FILE* stream is provided then table will be output to this. For
/// example it could be stdout or stderr.
	if (tableconfig->stream) {
		stream = tableconfig->stream;

/// If a neither a filename, name, or stream is defined in the
/// TABLECONFIG then the output is written to a filename which is the
/// name in the PMX object extended with a count of the tables written
/// so far.
	} else {
		if (!filename) {
			*tablecount += 1;
			filename = pmx_filename;
			sprintf(ext, "table.%i" OPENPMX_TABLEFILE, *tablecount);
		}

		char fname[PATH_MAX + NAME_MAX] = "";
		if (filename && filename[0]) { /* if its either 0 or "" we dont use it */
			strcpy(fname, filename);
			strcat(fname, ".");
		}
		strcat(fname, ext);
		stream = fopen(fname, "w");
	}
	assert(stream);

	var ret = (TABLE) {
		.idata = idata,
		.advanfuncs = advanfuncs,
		.popmodel = *popmodel,
		.tableconfig = tableconfig,

		.zeroerr = { },

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

/// The table fields are indicated by a string which will be tokenized
/// with spaces, tabs or comma.
	/* parse fields string and write the header */
	/* we have to tokenize a copy of the fields */
	char* token;
	char* rest = ret.fields;
	while ((token = strtok_r(rest, ", \t", &rest))) {
		strip_firstlast_space(token);
		vector_append(ret.fieldnames, token);
		fprintf(stream, OPENPMX_HEADER_FORMAT, token);
	}
	fputc('\n', stream);

	return ret;
}

static void table_close(TABLE* const table)
{
	/* close the stream if we opened it */
	if (table->stream && !table->tableconfig->stream)
		fclose(table->stream);
	table->stream = 0;

	/* cleanup */
	vector_free(table->fieldnames);
	free(table->fields);
}

/// Table generation is single threaded because no advancing takes 
/// place, only the saved state is used at each step. So it should be I/O
/// limited and not cpu limited. It would be very hard to make it
/// multithreaded because the rows have to be output in the correct
/// order.
static int table_row(TABLE* const table)
{
	assert(table);
	let idata = table->idata;
	let advanfuncs = table->advanfuncs;

	/* for first iteration setup, everything is 0 */
	if (table->record == 0) {
		table->individ = 0;
		table->individ_i = 0;

	/* we are already handling a normal row */
	} else {
		bool firstonly = false;
		if (table->tableconfig)
			firstonly = table->tableconfig->firstonly;
		let nrecord = idata->individ[table->individ].nrecord;

		/// If firstonly flag in the TABLECONFIG is set then only the
		/// first record of each individual in output to the table.
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
			return 0;
		}
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

	/* Table entries for each RECORD */
	if (strcmp(name, "yhat") == 0 ||
		strcmp(name, "YHAT") == 0)
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

	/* Table entries for THETA, ETA, STATE, or ERR values */
	const int off = _offset1 ? 1 : 0;
	var i = -1;
	if (sscanf(name, "theta[%i]", &i) == 1 ||
		sscanf(name, "THETA[%i]", &i) == 1 ||
		sscanf(name, "theta%i", &i) == 1 ||
		sscanf(name, "THETA%i", &i) == 1) {
		i -= off;
		if (i >= 0 && i < idata->ntheta)
			return table->popmodel.theta[i];
		if (table->individ == 0 && table->individ_i == 0)
			warning(0, "Accessing (\"%s\") outside ntheta (%i)\n", name, idata->ntheta + off);
		return NAN;
	}
	if (sscanf(name, "eta[%i]", &i) == 1 ||
		sscanf(name, "ETA[%i]", &i) == 1 ||
		sscanf(name, "eta%i", &i) == 1 ||
		sscanf(name, "ETA%i", &i) == 1) {
		i -= off;
		if (i >= 0 && i < idata->nomega)
			return table->eta[i];
		if (table->individ == 0 && table->individ_i == 0)
			warning(0, "Accessing (\"%s\") outside nomega (%i)\n", name, idata->nomega + off);
		return NAN;
	}
	if (sscanf(name, "a[%i]", &i) == 1 ||
		sscanf(name, "state[%i]", &i) == 1 ||
		sscanf(name, "A[%i]", &i) == 1 || 
		sscanf(name, "STATE[%i]", &i) == 1) {
		i -= off;
		if (i >= 0 && i < idata->nstate)
			return table->state[i];	
		if (table->individ == 0 && table->individ_i == 0)
			warning(0, "Accessing (\"%s\") outside nstate (%i)\n", name, idata->nstate + off);
		return NAN;
	}
	if (sscanf(name, "err[%i]", &i) == 1) {
		i -= off;
		if (i >= 0 && i < idata->nsigma)
			return table->err[i];
		if (table->individ == 0 && table->individ_i == 0)
			warning(0, "Accessing (\"%s\") outside nsigma (%i)\n", name, idata->nsigma + off);
		return NAN;
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
			let name = table.fieldnames.ptr[i];
			let val = table_value(&table, name, pmx->data._offset1);
			fprintf(f, OPENPMX_TABLE_FORMAT, val);
		}
		fputc('\n', f);
	}
	table_close(&table);
}

