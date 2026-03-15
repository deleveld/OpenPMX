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

#include <math.h>

#include "openpmx.h"
#include "defines.h"
#include "popmodel.h"
#include "utils/c22.h"
#include "utils/vector.h"
#include "utils/getdelim.h"
#include "utils/errctx.h"

typedef VECTOR(const char*) TOKENVEC;

typedef struct {
	enum {
		COL_IGNORE = 0,
		COL_THETA,
		COL_OMEGA,
		COL_SIGMA
	} kind;
	double value;
	int fixed;
	
	/* for theta */
    int num;        /* THETA or SIGMA number (1-based) */
	double theta_upper;
	double theta_lower;

	/* for omega */
    int row, col;   /* OMEGA indices (1-based) */
	int blocknum;
} COLINFO;

typedef VECTOR(COLINFO) EXTCOLS;

static COLINFO parse_col_name(const char *name, ERRCTX* errctx)
{
	COLINFO ci = { 0 };

	if (strncmp(name, "THETA", 5) == 0) {
		ci.kind = COL_THETA;
		const char *p = name + 5;
		char *end;
		long num = strtol(p, &end, 10);
		if (end == p || *end != '\0') {
			add_errctx(errctx, "%s: cannot parse theta index in \"%s\"", __func__, name);
			return ci;
		}
		if (num <= 0) {
			add_errctx(errctx, "%s: theta index must be >= 1 in \"%s\"", __func__, name);
			return ci;
		}
		if (num >= OPENPMX_THETA_MAX) {
			add_errctx(errctx, "%s: theta index %ld exceeds max (%d) in \"%s\"", __func__, num, OPENPMX_THETA_MAX, name);
			return ci;
		}
		ci.num = (int)num;

	} else if (strncmp(name, "SIGMA(", 6) == 0) {
		ci.kind = COL_SIGMA;
		const char *p = name + 6;
		char *end;
		long row = strtol(p, &end, 10);
		if (end == p || *end != ',') {
			add_errctx(errctx, "%s: cannot parse sigma indices in \"%s\"", __func__, name);
			return ci;
		}
		p = end + 1;
		long col = strtol(p, &end, 10);
		if (end == p || (*end != ')' && *end != '\0')) {
			add_errctx(errctx, "%s: cannot parse sigma col index in \"%s\"", __func__, name);
			return ci;
		}
		if (row <= 0 || col <= 0) {
			add_errctx(errctx, "%s: sigma indices must be >= 1 in \"%s\"", __func__, name);
			return ci;
		}
		if (row >= OPENPMX_SIGMA_MAX || col >= OPENPMX_SIGMA_MAX) {
			add_errctx(errctx, "%s: sigma index exceeds max (%d) in \"%s\"", __func__, OPENPMX_SIGMA_MAX, name);
			return ci;
		}
		if (row != col) {
			add_errctx(errctx, "%s: sigma not on diagonal in \"%s\"", __func__, name);
			return ci;
		}
		ci.row = (int)row;
		ci.col = (int)col;
		ci.num = ci.row;

	} else if (strncmp(name, "OMEGA(", 6) == 0) {
		ci.kind = COL_OMEGA;
		const char *p = name + 6;
		char *end;
		long row = strtol(p, &end, 10);
		if (end == p || *end != ',') {
			add_errctx(errctx, "%s: cannot parse omega indices in \"%s\"\n", __func__, name);
			return ci;
		}
		p = end + 1;
		long col = strtol(p, &end, 10);
		if (end == p || (*end != ')' && *end != '\0')) {
			add_errctx(errctx, "%s: cannot parse omega col index in \"%s\"\n", __func__, name);
			return ci;
		}
		if (row <= 0 || col <= 0) {
			add_errctx(errctx, "%s: omega indices must be >= 1 in \"%s\"\n", __func__, name);
			return ci;
		}
		if (row >= OPENPMX_OMEGA_MAX) {
			add_errctx(errctx, "%s: omega row %ld exceeds max (%d) in \"%s\"\n", __func__, row, OPENPMX_OMEGA_MAX, name);
			return ci;
		}
		if (col >= OPENPMX_OMEGA_MAX) {
			add_errctx(errctx, "%s: omega col %ld exceeds max (%d) in \"%s\"\n", __func__, col, OPENPMX_OMEGA_MAX, name);
			return ci;
		}
		if (col > row) {
			add_errctx(errctx, "%s: omega not lower triangular in \"%s\"\n", __func__, name);
			return ci;
		}
		ci.row = (int)row;
		ci.col = (int)col;
	}
	return ci;
}

/* modifies line and makes a vector of pointers to tokens */
static void get_delim_token_ptrs(char* line, TOKENVEC* namevec)
{
	char *saveptr;
	let delim = " \t\r\n";
	char *token = strtok_r(line, delim, &saveptr); 
	while (token != NULL) {
		vector_append(*namevec, token);
		token = strtok_r(NULL, delim, &saveptr);
	}
}

static void parse_extcols_line(char *line,
							   EXTCOLS* extcols,
							   ERRCTX* errctx,
							   const char* filename,
							   const int linenum)
{
	/* read in a line of values */
	TOKENVEC vals = { 0 };
	get_delim_token_ptrs(line, &vals);

	/* error if number of values detected didnt match */
	if (vals.size != extcols->size) {
		add_errctx(errctx, "%s: ext file \"%s\" line %i number fields (%i) does not match header (%i)\n",
				   __func__, filename, linenum, vals.size, extcols->size);

	/* copy over values if no errors */
	} else {

		/* basic value */
		let iteration = atol(vals.ptr[0]);
		if (iteration > 0 || iteration == -1000000000) {
			forvector(i, vals) 
				extcols->rawptr[i].value = atof(vals.ptr[i]);

		/* lower limit */
		} else if (iteration == OPENPMX_EXTFILE_LOWER_BOUNDS) {
			forvector(i, vals) 
				extcols->rawptr[i].theta_lower = atof(vals.ptr[i]);

		/* upper limit */
		} else if (iteration == OPENPMX_EXTFILE_UPPER_BOUNDS) {
			forvector(i, vals) 
				extcols->rawptr[i].theta_upper = atof(vals.ptr[i]);

		/* fixed */
		} else if (iteration == OPENPMX_EXTFILE_EXTRA_FIXED) {
			forvector(i, vals) 
				extcols->rawptr[i].fixed = atoi(vals.ptr[i]);

		/* blocknum */
		} else if (iteration == OPENPMX_EXTFILE_OMEGA_BLOCKS) {
			forvector(i, vals) 
				extcols->rawptr[i].blocknum = atoi(vals.ptr[i]);
		}
	}
}

static void read_extcols(const char *filename, EXTCOLS* extcols, ERRCTX* errctx)
{
	char *line = NULL;
	char* header_text = 0;
	TOKENVEC header_elem = { 0 };

    FILE *stream = fopen(filename, "r");
    if (!stream) {
		add_errctx(errctx, "%s: cannot open \"%s\"\n", __func__, filename);
		goto failed;
	}

	/* read in header */
    size_t capacity = 0;
    var len = openpmx_getdelim(&line, &capacity, '\n', stream);
    if (len == -1) {
		add_errctx(errctx, "%s: ext file \"%s\" cannot read header\n", __func__, filename);
		goto failed;
	}
	var linenum = 1;

	/* read and translate header */
	header_text = strdup(line);
	get_delim_token_ptrs(header_text, &header_elem);
	forvector(i, header_elem) {
		let name = header_elem.ptr[i];
		let colinfo = parse_col_name(name, errctx);
		if (errctx->len) {
			add_errctx(errctx, "%s: ext file \"%s\" header parse fail\n", __func__, filename);
			goto failed;
		}
		vector_append(*extcols, colinfo);
	}

	/* read in lines, parse them, check for errors */
    while ((len = openpmx_getdelim(&line, &capacity, '\n', stream)) != -1) {
		++linenum;
		parse_extcols_line(line, extcols, errctx, filename, linenum);
		if (errctx->len)
			break;
	}

failed:
	/* cleanup stream, line and header */
	if (stream)
		fclose(stream);
	free(line);
	vector_free(header_elem);
	free(header_text);
}

static OPENPMX popparams_init_from_file(const char* filename, ERRCTX* errctx)
{
    OPENPMX pmx;
    memset(&pmx, 0, sizeof(pmx));
    forcount(i, OPENPMX_THETA_MAX) 
		pmx.theta[i].type = THETA_INVALID;
    forcount(i, OPENPMX_SIGMA_MAX)
		pmx.sigma[i] = 0;
    forcount(i, OPENPMX_OMEGABLOCK_MAX)
		pmx.omega[i].type = OMEGA_INVALID;

	EXTCOLS extcols = { 0 };
	read_extcols(filename, &extcols, errctx);

	/* fail with any lower level error */
	if (errctx->len)
		goto failed;

	double omega[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX] = { 0 };
	int omegafixed[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX] = { 0 };
	int omegablocknum[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX] = { 0 };

	forvector(i, extcols) {
		let v = &extcols.ptr[i];

		if (v->kind == COL_THETA) {
			let thetaidx = v->num - 1;
			pmx.theta[thetaidx].value = v->value;
			pmx.theta[thetaidx].lower = v->theta_lower;
			pmx.theta[thetaidx].upper = v->theta_upper;
			pmx.theta[thetaidx].type = v->fixed ? FIXED : ESTIMATE;

		} else if (v->kind == COL_SIGMA) {
			let sigmaidx = v->num - 1;
			var value = v->value;
			if (v->fixed)
				value = -fabs(value);
			pmx.sigma[sigmaidx] = value;

		} else if (v->kind == COL_OMEGA) {
			let rowidx = v->row - 1;
			let colidx = v->col - 1;
			omega[rowidx][colidx] = omega[colidx][rowidx] = v->value;
			omegafixed[rowidx][colidx] = omegafixed[colidx][rowidx] = v->fixed;
			omegablocknum[rowidx][colidx] = omegablocknum[colidx][rowidx] = v->blocknum;
		}
	}
	/* write the pmx omega block structure */
	/* Find nomega: largest row or col index seen in any OMEGA column.
	 * Indices are 1-based in the ext file.                            */
	int nomega = 0;
	forvector(i, extcols) {
		let v = &extcols.ptr[i];
		if (v->kind == COL_OMEGA) {
			if (v->row > nomega) nomega = v->row;
			if (v->col > nomega) nomega = v->col;
		}
	}

	/* Walk the diagonal (0-based) and emit one pmx.omega[] block per
	 * run of consecutive diagonal entries sharing the same block number. */
	int bk = 0;   /* index into pmx.omega[]  */
	int d  = 0;   /* 0-based diagonal offset */

	while (d < nomega && bk < OPENPMX_OMEGABLOCK_MAX) {
		int bn  = omegablocknum[d][d];

		/* How many consecutive diagonal entries share this block number? */
		int dim = 0;
		while (d + dim < nomega && omegablocknum[d + dim][d + dim] == bn)
			dim++;

		/* Determine block type:
		 *   OMEGA_SAME  – every diagonal element has fixed == 2
		 *   OMEGA_BLOCK – at least one off-diagonal has a non-zero blocknum
		 *   OMEGA_DIAG  – otherwise                                       */
		int all_same    = 1;
		int has_offdiag = 0;
		for (int r = 0; r < dim; r++) {
			if (omegafixed[d + r][d + r] != 2)
				all_same = 0;
			for (int c = 0; c < r; c++)
				if (omegablocknum[d + r][d + c] != 0)
					has_offdiag = 1;
		}
		int type;
		if (all_same && dim > 0)
			type = OMEGA_SAME;
		else if (has_offdiag)
			type = OMEGA_BLOCK;
		else
			type = OMEGA_DIAG;

		pmx.omega[bk].type = type;
		pmx.omega[bk].ndim = dim;

		/* Pack values[] the way popmodel_init() expects to read them:
		 *
		 *   OMEGA_DIAG:  one value per diagonal element.
		 *                Negative value encodes a fixed diagonal.
		 *
		 *   OMEGA_BLOCK: lower-triangle row-major per row —
		 *                off-diagonal elements first, then the diagonal.
		 *                Negative diagonal encodes a fixed diagonal.
		 *
		 *   OMEGA_SAME:  no values needed; popmodel_init() copies from
		 *                the immediately preceding block.               */
		if (type == OMEGA_DIAG) {
			int n = 0;
			for (int r = 0; r < dim; r++) {
				double val = omega[d + r][d + r];
				if (omegafixed[d + r][d + r] == 1)
					val = -fabs(val);
				pmx.omega[bk].values[n++] = val;
			}

		} else if (type == OMEGA_BLOCK) {
			int n = 0;
			for (int r = 0; r < dim; r++) {
				/* off-diagonal entries for this row */
				for (int c = 0; c < r; c++)
					pmx.omega[bk].values[n++] = omega[d + r][d + c];
				/* diagonal entry — negative if fixed */
				double val = omega[d + r][d + r];
				if (omegafixed[d + r][d + r] == 1)
					val = -fabs(val);
				pmx.omega[bk].values[n++] = val;
			}

		} /* OMEGA_SAME: leave values[] zeroed, popmodel_init copies */

		d  += dim;
		bk += 1;
	}
	
failed:
    vector_free(extcols);
    return pmx;
}

int pmx_reload_popparam(OPENPMX* dest, RELOADCONFIG* args)
{
	ERRCTX errctx = { 0 };

	/* if no filename passed in, make one from the destitation name */
	const char* filename = args->filename;
	char* filename_buffer = 0;
	if (!filename) {
		if (!dest->filename) {
			add_errctx(&errctx, "%s: no default ext filename\n", __func__);
			goto failed;
		}
		/* construct the default filename, where we expect the ext file
		 * would normally be */
		let len = strlen(dest->filename) + strlen(OPENPMX_EXTFILE);
		filename_buffer = mallocvar(char, len + 1);
		strcpy(filename_buffer, dest->filename);
		strcat(filename_buffer, OPENPMX_EXTFILE);
		filename = filename_buffer;
	}
	/* read in the source from the ext file */
    var src = popparams_init_from_file(filename, &errctx);
	if (errctx.len)
		goto failed;

	/* we have to check whether we do the copy in its entirety before
	 * even starting, so we dont copy over a half initialized object */
	if (!args->force) {
		forcount(i, OPENPMX_THETA_MAX) {
			if (src.theta[i].type == THETA_INVALID)
				break;
			if (dest->theta[i].lower != src.theta[i].lower ||
				dest->theta[i].upper != src.theta[i].upper) {
				add_errctx(&errctx, "%s: theta bounds mismatch\n", __func__);
				goto failed;
			}
			if (dest->theta[i].type != src.theta[i].type) {
				add_errctx(&errctx, "%s: theta type mismatch\n", __func__);
				goto failed;
			}
			if (src.theta[i].value <= dest->theta[i].lower ||
				src.theta[i].value >= dest->theta[i].upper) {
				add_errctx(&errctx, "%s: theta value outside range\n", __func__);
				goto failed;
			}
		}

		/* sigma cant fail */

		forcount(i, OPENPMX_OMEGABLOCK_MAX) {
			if (src.omega[i].type == OMEGA_INVALID)
				break;
			if (dest->omega[i].type != src.omega[i].type) {
				add_errctx(&errctx, "%s: omega block type mismatch\n", __func__);
				goto failed;
			}
			if (dest->omega[i].ndim != src.omega[i].ndim) {
				add_errctx(&errctx, "%s: omega block ndim mismatch\n", __func__);
				goto failed;
			}
		}
	}

	/* success path */
	/* actually copy over the population parameters */
	int i;
	for (i=0; i<OPENPMX_THETA_MAX; i++) {
		if (src.theta[i].type == THETA_INVALID)
			break;
		dest->theta[i] = src.theta[i];
	}
	let ntheta = i;

	/* last non-zero sigma is the last one */
	var nsigma = 0;
	for (i=0; i<OPENPMX_SIGMA_MAX; i++) {
		if (src.sigma[i] != 0.)
			nsigma = i + 1;
	}
	for (i=0; i<nsigma; i++)
		dest->sigma[i] = src.sigma[i];
	
	for (i=0; i<OPENPMX_OMEGABLOCK_MAX; i++) {
		if (src.omega[i].type == OMEGA_INVALID)
			break;
		dest->omega[i] = src.omega[i];
	}
	let nblock = i;

	/* if we dont preserve old values then overwrite other paramaters */
	if (!args->preserve) {
		for (int i=ntheta; i<OPENPMX_THETA_MAX; i++)
			memset(&dest->theta[i], 0, sizeof(dest->theta[i]));

		for (int i=nsigma; i<OPENPMX_SIGMA_MAX; i++)
			dest->sigma[i] = 0.;

		for (int i=nblock; i<OPENPMX_OMEGABLOCK_MAX; i++)
			memset(&dest->omega[i], 0, sizeof(dest->omega[i]));
	}

	/* test by making a POPMODEL which catches more errors */
	let popmodel = popmodel_init(dest, &errctx);
	if (errctx.len)
		goto failed;
	if (!args->silent) {
		fprintf(stdout, "popparam loaded \"%s\"\n", filename);
		popmodel_information(0, &popmodel, 0);
	}

	/* cleanup */
	/* return with success */
	dest->result = (PMXRESULT) { 0 };
	if (filename_buffer)
		free(filename_buffer);
	return 0;

	/* failure path */
failed:
	free(filename_buffer);

	if (errctx.len)
		fprintf(stderr, "%s", errctx.errmsg);

	/* fail if we are not optional */
	if (!args->optional)
		exit(EXIT_FAILURE);
		
	return 1;
}


