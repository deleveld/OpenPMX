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
#include "utils/getdelim.h"

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
	const double* err;
} TABLE;

#ifndef NAME_MAX
#define NAME_MAX 256
#endif

static void table_close(TABLE* const table)
{
	assert(table);
	
	/* close the stream if there is one and if we opened it */
	if (table->stream &&
		table->stream != stdout &&
		table->stream != stderr &&
		(!table->tableconfig || !table->tableconfig->stream))
		fclose(table->stream);
	table->stream = 0;

	/* cleanup */
	vector_free(table->fieldnames);
	free(table->fields);
}

static bool check_for_whitespace(const char* p)
{
	let len = strlen(p);
	forcount(i, len) {
		if (isspace(p[i]))
			return true;
	}
	return false;
}

static TABLE table_open(const IDATA* const idata,
						const ADVANFUNCS* const advanfuncs,
						const POPMODEL* const popmodel,
						const char* pmx_filename,
						int* tablecount,
						const char* _fields,
						const TABLECONFIG* const tableconfig,
						ERRCTX* errctx)
{
	FILE* stream = 0;
	
	if (tableconfig) {
/// If a .filename is defined in the TABLECONFIG then the output is
/// written to that filename.  
		/* write to user define filename */
		if (tableconfig->filename) {
			if (tableconfig->name) {
				errctx_add(errctx, "%s: table has both filename \"%s\" and name \"%s\"\n",
						   __func__, tableconfig->filename, tableconfig->name);
				goto done_file_open;
			}
			if (tableconfig->stream) {
				errctx_add(errctx, "%s: table has both filename \"%s\" and stream\n",
						   __func__, tableconfig->filename);
				goto done_file_open;
			}
			stream = fopen(tableconfig->filename, "w");
			if (!stream) {
				errctx_add(errctx, "%s: open table file failed \"%s\"\n", __func__, tableconfig->filename);
				goto done_file_open;
			}
/// If a .name is defined in the TABLECONFIG then the output is
/// written to a filename which is the name in the PMX object extended
/// with the given name.   
		} else if (tableconfig->name) {
			char fname[PATH_MAX + NAME_MAX];
			if (tableconfig->stream) {
				errctx_add(errctx, "%s: table has both name \"%s\" and stream\n",
						   __func__, tableconfig->name);
				goto done_file_open;
			}
			snprintf(fname, sizeof(fname), "%s%s%s",
				pmx_filename ? pmx_filename : "",
				pmx_filename ? "." : "",
				tableconfig->name);
			stream = fopen(fname, "w");
			if (!stream) {
				errctx_add(errctx, "%s: table open failed \"%s\"\n",
						   __func__, fname);
				goto done_file_open;
			}
/// if a .stream (FILE*) is provided then table will be output to it.
/// For example it could be stdout or stderr. If given in this way then
/// fclose() is not called on the stream when the table is closed.
		} else if (tableconfig->stream) {
			stream = tableconfig->stream;
		}
	}
/// If no TABLECONFIG, then see if we have a filename in PMX then make
/// extension with numbered tables.
	if (!stream) { 
		if (pmx_filename) {
			char fname[PATH_MAX + NAME_MAX];
			*tablecount += 1;
			snprintf(fname, sizeof(fname), "%s.table.%i%s",
				pmx_filename,
				*tablecount,
				OPENPMX_TABLEFILE);
			stream = fopen(fname, "w");
			if (!stream) {
				errctx_add(errctx, "%s: open table file failed \"%s\"\n", __func__, fname);
				goto done_file_open;
			}

/// If no TABLECONFIG and no PMX filename, then use stdout
		} else
			stream = stdout;
	}

done_file_open:
	if (!stream) {
		errctx_add(errctx, "%s: no stream\n", __func__);
		return (TABLE) { 0 };
	}
	
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
		.err = 0,
	};

/// The table fields are indicated by a string which will be tokenized
/// with whitespace or comma.
	VECPTR tokens = { 0 };
	let err = get_delim_tokens(ret.fields, &tokens, GET_DELIM_SEP_ANY);
	if (err) {
		errctx_add(errctx, "%s: table names tokens failed\n", __func__);
		table_close(&ret);
		return (TABLE) { 0 };
	}
	
	/* check for invalid table fields */
	forvector_val(p, tokens) {
		if (check_for_whitespace(p)) {
			errctx_add(errctx, "%s: table field \"%s\" is invalid\n", __func__, p);
			table_close(&ret);
			return (TABLE) { 0 };
		}
	}

	/* table fields are OK, use them */
	vector_appendn(ret.fieldnames, tokens.ptr, tokens.size);
	vector_free(tokens);

	/* write the header out */
	forvector_val(fieldname, ret.fieldnames) 
		fprintf(stream, OPENPMX_HEADER_FORMAT, fieldname);
	fputc('\n', stream);

	return ret;
}

/// Table generation is single threaded because no advancing takes 
/// place, only the saved state is used at each step. Thus it is I/O
/// limited and not cpu limited.
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
		var firstonly = false;
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

/// All of the fields of the data record are accessible.
	let recordinfo = &advanfuncs->recordinfo;
	let recordfields = &recordinfo->dataconfig->recordfields;
	var recordoffset = structinfo_find_offset(name, recordfields);
	if (recordoffset >= 0)
		return DATA_FIELD(table->record, recordoffset);

/// All of the fields of $IMODEL(...) are accessible.
	let modelfields = &advanfuncs->advanconfig->imodelfields;
	let modeloffset = structinfo_find_offset(name, modelfields);
	if (modeloffset >= 0)
		return DATA_FIELD(table->imodel, modeloffset);

/// All of the fields of $PREDICT(...) are accessible.
	let predictfields = &advanfuncs->advanconfig->predictfields;
	let predictvarsoffset = structinfo_find_offset(name, predictfields);
	if (predictvarsoffset >= 0)
		return DATA_FIELD(table->predictvars, predictvarsoffset);

/// The follwing fields are accesible: YHAT, YHATVAR, PRED, OBJ, INEVAL,
/// EVID, MDV, EVID, AMT, RATE, and CMT.
	/* Table entries for each RECORD */
	if (strcmp(name, "YHAT") == 0)
		return individ->yhat[table->individ_i];
	if (strcmp(name, "YHATVAR") == 0)
		return individ->yhatvar[table->individ_i];
	if (strcmp(name, "PRED") == 0)
		return individ->pred[table->individ_i];
	if (strcmp(name, "OBJ") == 0)
		return individ->iobjfn;
	if (strcmp(name, "INEVAL") == 0)
		return individ->ineval; 
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
	if (sscanf(name, "THETA(%i)", &i) == 1) {
		i -= off;
		if (i >= 0 && i < idata->ntheta)
			return table->popmodel.theta[i];
		if (table->individ == 0 && table->individ_i == 0)
			warning(0, "access \"%s\" outside ntheta (%i)\n", name, idata->ntheta + off);
		return NAN;
	}
	if (sscanf(name, "ETA(%i)", &i) == 1) {
		i -= off;
		let eta = individ->eta;
		if (i >= 0 && i < idata->nomega)
			return eta[i];
		if (table->individ == 0 && table->individ_i == 0)
			warning(0, "access \"%s\" outside nomega (%i)\n", name, idata->nomega + off);
		return NAN;
	}
	if (sscanf(name, "A(%i)", &i) == 1) {
		i -= off;
		let state = individ->istate + table->individ_i * idata->nstate;
		if (i >= 0 && i < idata->nstate)
			return state[i];	
		if (table->individ == 0 && table->individ_i == 0)
			warning(0, "access \"%s\" outside nstate (%i)\n", name, idata->nstate + off);
		return NAN;
	}
	if (sscanf(name, "ERR(%i)", &i) == 1) {
		i -= off;
		if (i >= 0 && i < idata->nsigma)
			return table->err[i];
		if (table->individ == 0 && table->individ_i == 0)
			warning(0, "access \"%s\" outside nsigma (%i)\n", name, idata->nsigma + off);
		return NAN;
	}
	
	warning(0, "Cannot find table field \"%s\", using NAN.\n", name);
	return NAN;
}

void pmx_table(OPENPMX* pmx,
			   const char* fields,
			   const TABLECONFIG* const tableconfig)
{
	pmxstate_ensure(pmx);

	ERRCTX errctx = { 0 };
	let popmodel = popmodel_init(pmx->theta, pmx->omega, pmx->sigma, &errctx);
	if (errctx.len)
		fatal(0, "%s: %s", __func__, errctx.errmsg);

	TABLE table = table_open(&pmx->state->idata,
							 pmx->state->advanfuncs,
							 &popmodel,
							 pmx->filename,
							 &pmx->state->tablecount,
							 fields,
							 tableconfig,
							 &errctx);
	if (errctx.len)
		fatal(0, "%s: %s", __func__, errctx.errmsg);

	FILE* f = table.stream;
	assert(f);

	/* write out fields */
	while (table_row(&table)) {
		forvector_val(name, table.fieldnames) {
			let val = table_value(&table, name, pmx->data._offset1);
			fprintf(f, OPENPMX_TABLE_FORMAT, val);
		}
		fputc('\n', f);
	}
	table_close(&table);
}

