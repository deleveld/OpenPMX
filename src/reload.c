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
#include "omegafixed.h"
#include "idata.h"
#include "print.h"
#include "pmxstate.h"
#include "advan/advan.h"
#include "utils/c22.h"
#include "utils/vector.h"
#include "utils/getdelim.h"
#include "utils/errctx.h"

/// This file implements a function that modifies an OPENPMX object
/// from information from a .ext file. It is available in openpmxtran
/// as `reload()`. This loads the population paramaters, setting their
/// structure to that in the file.

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

static COLINFO parse_col_name(const char *name, ERRCTX* errctx)
{
	COLINFO ci = { 0 };

	if (strncmp(name, "THETA", 5) == 0) {
		ci.kind = COL_THETA;
		const char *p = name + 5;
		char *end;
		long num = strtol(p, &end, 10);
		if (end == p || *end != '\0') {
			errctx_add(errctx, "%s: cannot parse theta index in \"%s\"\n", __func__, name);
			return ci;
		}
		if (num <= 0) {
			errctx_add(errctx, "%s: theta index must be >= 1 in \"%s\"\n", __func__, name);
			return ci;
		}
		if (num >= OPENPMX_THETA_MAX) { /* num is 1-based, conversion to 0-based occurs later */
			errctx_add(errctx, "%s: theta index %ld exceeds max (%d) in \"%s\"\n", __func__, num, OPENPMX_THETA_MAX, name);
			return ci;
		}
		ci.num = (int)num;

	} else if (strncmp(name, "SIGMA(", 6) == 0) {
		ci.kind = COL_SIGMA;
		const char *p = name + 6;
		char *end;
		long row = strtol(p, &end, 10);
		if (end == p || *end != ',') {
			errctx_add(errctx, "%s: cannot parse sigma indices in \"%s\"\n", __func__, name);
			return ci;
		}
		p = end + 1;
		long col = strtol(p, &end, 10);
		if (end == p || (*end != ')' && *end != '\0')) {
			errctx_add(errctx, "%s: cannot parse sigma col index in \"%s\"\n", __func__, name);
			return ci;
		}
		if (row <= 0 || col <= 0) {
			errctx_add(errctx, "%s: sigma indices must be >= 1 in \"%s\"\n", __func__, name);
			return ci;
		}
		if (row >= OPENPMX_SIGMA_MAX || col >= OPENPMX_SIGMA_MAX) { /* num is 1-based, conversion to 0-based occurs later */
			errctx_add(errctx, "%s: sigma index exceeds max (%d) in \"%s\"\n", __func__, OPENPMX_SIGMA_MAX, name);
			return ci;
		}
		if (row != col) {
			errctx_add(errctx, "%s: sigma not on diagonal in \"%s\"\n", __func__, name);
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
			errctx_add(errctx, "%s: cannot parse omega indices in \"%s\"\n", __func__, name);
			return ci;
		}
		p = end + 1;
		long col = strtol(p, &end, 10);
		if (end == p || (*end != ')' && *end != '\0')) {
			errctx_add(errctx, "%s: cannot parse omega col index in \"%s\"\n", __func__, name);
			return ci;
		}
		if (row <= 0 || col <= 0) {
			errctx_add(errctx, "%s: omega indices must be >= 1 in \"%s\"\n", __func__, name);
			return ci;
		}
		if (row >= OPENPMX_OMEGA_MAX) { /* num is 1-based, conversion to 0-based occurs later */
			errctx_add(errctx, "%s: omega row %ld exceeds max (%d) in \"%s\"\n", __func__, row, OPENPMX_OMEGA_MAX, name);
			return ci;
		}
		if (col >= OPENPMX_OMEGA_MAX) { /* num is 1-based, conversion to 0-based occurs later */
			errctx_add(errctx, "%s: omega col %ld exceeds max (%d) in \"%s\"\n", __func__, col, OPENPMX_OMEGA_MAX, name);
			return ci;
		}
		if (col > row) {
			errctx_add(errctx, "%s: omega not lower triangular in \"%s\"\n", __func__, name);
			return ci;
		}
		ci.row = (int)row;
		ci.col = (int)col;
	}
	return ci;
}

/* modifies line and makes a vector of pointers to tokens */
static void get_delim_token_ptrs(char* line,
								 TOKENVEC* namevec)
{
	char *saveptr;
	let delim = " \t\r\n";
	char *token = strtok_r(line, delim, &saveptr); 
	while (token != NULL) {
		vector_append(*namevec, token);
		token = strtok_r(NULL, delim, &saveptr);
	}
}

typedef VECTOR(COLINFO) EXTCOLS;

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
		errctx_add(errctx, "%s: ext file \"%s\" line %i number fields (%i) does not match header (%i)\n",
				   __func__, filename, linenum, vals.size, extcols->size);

	/* copy over values if no errors */
	} else {

/// Each line of the .ext file is read and saved. This allows
/// continuation of terminated estimation runs from the population
/// paramaters from the last full iteration.
		/* basic value */
		let iteration = atol(vals.ptr[0]);
		if (iteration > 0 || iteration == -1000000000) {
			forvector(i, vals) 
				extcols->mutptr[i].value = atof(vals.ptr[i]);

		/* lower limit */
		} else if (iteration == OPENPMX_EXTFILE_LOWER_BOUNDS) {
			forvector(i, vals) 
				extcols->mutptr[i].theta_lower = atof(vals.ptr[i]);

		/* upper limit */
		} else if (iteration == OPENPMX_EXTFILE_UPPER_BOUNDS) {
			forvector(i, vals) 
				extcols->mutptr[i].theta_upper = atof(vals.ptr[i]);

		/* fixed */
		} else if (iteration == OPENPMX_EXTFILE_EXTRA_FIXED) {
			forvector(i, vals) 
				extcols->mutptr[i].fixed = atoi(vals.ptr[i]);

		/* blocknum */
		} else if (iteration == OPENPMX_EXTFILE_OMEGA_BLOCKS) {
			forvector(i, vals) 
				extcols->mutptr[i].blocknum = atoi(vals.ptr[i]);
		}
	}
	vector_free(vals);
}
	
static void read_extcols(const char *filename, EXTCOLS* extcols, ERRCTX* errctx)
{
	char *line = NULL;
	char* header_text = 0;
	TOKENVEC header_elem = { 0 };

    FILE *stream = fopen(filename, "r");
    if (!stream) {
		errctx_add(errctx, "%s: cannot open \"%s\"\n", __func__, filename);
		goto failed;
	}

	/* read in header */
    size_t capacity = 0;
    var len = openpmx_getdelim(&line, &capacity, '\n', stream);
    if (len == -1) {
		errctx_add(errctx, "%s: ext file \"%s\" cannot read header\n", __func__, filename);
		goto failed;
	}
	var linenum = 1;

	/* read and translate header */
	header_text = strdup(line);
	get_delim_token_ptrs(header_text, &header_elem);
	forvector_val(name, header_elem) {
		let colinfo = parse_col_name(name, errctx);
		if (errctx->len) {
			errctx_add(errctx, "%s: ext file \"%s\" header parse fail\n", __func__, filename);
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

typedef struct {
	typeof(((OPENPMX){0}).theta) theta;
	typeof(((OPENPMX){0}).omega) omega;
	typeof(((OPENPMX){0}).sigma) sigma;
} RELOADPARAM;

static RELOADPARAM params_init_from_file(const char* filename,
										 ERRCTX* errctx)
{
    RELOADPARAM p = { 0 };
    forcount(i, OPENPMX_THETA_MAX) 
		p.theta[i].type = THETA_INVALID;
    forcount(i, OPENPMX_SIGMA_MAX)
		p.sigma[i] = 0;
    forcount(i, OPENPMX_OMEGABLOCK_MAX)
		p.omega[i].type = OMEGA_INVALID;

	/* read in coloumns and exit if error */
	EXTCOLS extcols = { 0 };
	read_extcols(filename, &extcols, errctx);
	if (errctx->len)
		goto failed;

	double omega[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX] = { 0 };
	OMEGAFIXED omegafixed[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX] = { 0 };
	int omegablocknum[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX] = { 0 };

	forvector_ptr(v, extcols) {
		if (v->kind == COL_THETA) {
			let thetaidx = v->num - 1;
			p.theta[thetaidx].value = v->value;
			p.theta[thetaidx].lower = v->theta_lower;
			p.theta[thetaidx].upper = v->theta_upper;
			p.theta[thetaidx].type = v->fixed ? FIXED : ESTIMATE;

		} else if (v->kind == COL_SIGMA) {
			let sigmaidx = v->num - 1;
			var value = v->value;
			if (v->fixed)
				value = -fabs(value);
			p.sigma[sigmaidx] = value;

		} else if (v->kind == COL_OMEGA) {
			let rowidx = v->row - 1;
			let colidx = v->col - 1;
			omega[rowidx][colidx] = omega[colidx][rowidx] = v->value;
			
			var fixval = omegafixed_decode(v->fixed);
			omegafixed[rowidx][colidx] = omegafixed[colidx][rowidx] = fixval; 
			omegablocknum[rowidx][colidx] = omegablocknum[colidx][rowidx] = v->blocknum;
		}
	}
	/* write the pmx omega block structure */
	/* Find nomega: largest row or col index seen in any OMEGA column.
	 * Indices are 1-based in the ext file.                            */
	int nomega = 0;
	forvector_ptr(v, extcols) {
		if (v->kind == COL_OMEGA) {
			if (v->row > nomega) nomega = v->row;
			if (v->col > nomega) nomega = v->col;
		}
	}

	/* Walk the diagonal (0-based) and emit one p.omega[] block per
	 * run of consecutive diagonal entries sharing the same block number. */
	var bk = 0;		/* index into p.omega[]  */
	var d = 0;		/* 0-based diagonal offset */

	while (d < nomega && bk < OPENPMX_OMEGABLOCK_MAX) {
		var bn = omegablocknum[d][d];

		/* How many consecutive diagonal entries share this block number? */
		var dim = 0;
		while (d + dim < nomega && omegablocknum[d + dim][d + dim] == bn)
			dim++;

		/* Determine block type:
		 *   OMEGA_SAME  – every diagonal element has fixed == 2
		 *   OMEGA_BLOCK – at least one off-diagonal has a non-zero blocknum
		 *   OMEGA_DIAG  – otherwise                                       */
		var all_same    = 1;
		var has_offdiag = 0;
		for (int r = 0; r < dim; r++) {
			if (omegafixed[d + r][d + r] != OMEGAFIXED_SAME)
				all_same = 0;
			for (int c = 0; c < r; c++)
				if (omegablocknum[d + r][d + c] != 0)
					has_offdiag = 1;
		}
		/* 0-dim block is error */
		if (dim == 0) {
			errctx_add(errctx, "%s: omega block size 0\n", __func__);
			goto failed;
		}
		
		
		var type = OMEGA_DIAG;
		if (all_same && dim > 0)
			type = OMEGA_SAME;
		else if (has_offdiag)
			type = OMEGA_BLOCK;

		p.omega[bk].type = type;
		p.omega[bk].ndim = dim;

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
				if (omegafixed[d + r][d + r] == OMEGAFIXED_FIXED)
					val = -fabs(val);
				p.omega[bk].values[n++] = val;
			}

		} else if (type == OMEGA_BLOCK) {
			int n = 0;
			for (int r = 0; r < dim; r++) {
				/* off-diagonal entries for this row */
				for (int c = 0; c < r; c++)
					p.omega[bk].values[n++] = omega[d + r][d + c];
				/* diagonal entry — negative if fixed */
				double val = omega[d + r][d + r];
				if (omegafixed[d + r][d + r] == OMEGAFIXED_FIXED)
					val = -fabs(val);
				p.omega[bk].values[n++] = val;
			}

		} /* OMEGA_SAME: leave values[] zeroed, popmodel_init copies */

		d  += dim;
		bk += 1;
	}
	
failed:
    vector_free(extcols);
    return p;
}

static int reload_popparam(OPENPMX* dest, RELOADCONFIG* args, POPMODEL* popmodel)
{
	ERRCTX errctx = { 0 };

	const char* filename = args->filename;
	char* filename_buffer = 0;
	if (!filename) {
		if (!dest->filename) {
			errctx_add(&errctx, "%s: no default ext filename\n", __func__);
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
    var src = params_init_from_file(filename, &errctx);
	if (errctx.len)
		goto failed;

	/// The default is that reloading of the population paramaters fails if: 
	///
	/// + Any theta bounds dont match
	/// + Any theta FIXED/ESTIMATE dont match
	/// + The structure of any the omega blocks differ
	/// + Any omega values fixed or estimated dont match
	/// + The number of sigma values differ
	/// + Any sigma values fixed or estimated dont match
	///
	/* we have to check whether we do the copy in its entirety before
	 * even starting, so we dont copy over a half initialized object */
	let s = popmodel_init(src.theta, src.omega, src.sigma, &errctx);
	if (errctx.len)
		goto failed;

	let d = popmodel_init(dest->theta, dest->omega, dest->sigma, &errctx);
	if (errctx.len)
		goto failed;

	if (!args->force) {
		var ntheta_bad = (s.ntheta != d.ntheta);
		var nblock_bad = (s.nblock != d.nblock);
		var nomega_bad = (s.nomega != d.nomega);
		var nsigma_bad = (s.nsigma != d.nsigma);
		if (ntheta_bad) {
			errctx_add(&errctx, "%s: theta count mismatch\n", __func__);
			goto failed;
		}
		if (nblock_bad) {
			errctx_add(&errctx, "%s: omega block count mismatch\n", __func__);
			goto failed;
		}
		if (nomega_bad) {
			errctx_add(&errctx, "%s: omega count mismatch\n", __func__);
			goto failed;
		}
		if (nsigma_bad) {
			errctx_add(&errctx, "%s: sigma count mismatch\n", __func__);
			goto failed;
		}

		forcount(i, s.ntheta) {
			if (d.lower[i] != s.lower[i] ||
				d.upper[i] != s.upper[i]) {
				errctx_add(&errctx, "%s: theta bounds mismatch\n", __func__);
				goto failed;
			}
			if (d.thetaestim[i] != s.thetaestim[i]) {
				errctx_add(&errctx, "%s: theta type mismatch\n", __func__);
				goto failed;
			}

			if (s.theta[i] <= d.lower[i] ||
				s.theta[i] >= d.upper[i]) {
				if (s.thetaestim[i] == ESTIMATE) {
					errctx_add(&errctx, "%s: estimated theta value outside range\n", __func__);
					goto failed;
				}
			}
		}

		forcount(i, s.nomega) {
			forcount(j, s.nomega) {
				if (s.omegafixed[i][j] !=  d.omegafixed[i][j]) {
					errctx_add(&errctx, "%s: omega fixed mismatch \n", __func__);
					goto failed;
				}
			}
		}

		forcount(i, s.nsigma) {
			if (s.sigmafixed[i] != d.sigmafixed[i]) {
				errctx_add(&errctx, "%s: sigma fixed mismatch\n", __func__);
				goto failed;
			}
		}
		
		forcount(i, s.nblock) {
			if (d.blocktype[i] != s.blocktype[i]) {
				errctx_add(&errctx, "%s: omega block type mismatch\n", __func__);
				goto failed;
			}
			if (d.blockdim[i] != s.blockdim[i]) {
				errctx_add(&errctx, "%s: omega block ndim mismatch\n", __func__);
				goto failed;
			}
		}
	}

	/* success path */ 
	/* actually copy over the population parameters, even the limits */
	memcpy(dest->theta, src.theta, s.ntheta*sizeof(THETATYPE));
	memcpy(dest->omega, src.omega, s.nblock*sizeof(OMEGABLOCKSTYPE));
	memcpy(dest->sigma, src.sigma, s.nsigma*sizeof(double));
	
	/* overwrite other paramaters */
	for (int i=s.ntheta; i<OPENPMX_THETA_MAX; i++)
		memset(&dest->theta[i], 0, sizeof(THETATYPE));

	for (int i=s.nsigma; i<OPENPMX_SIGMA_MAX; i++)
		dest->sigma[i] = 0.;

	for (int i=s.nblock; i<OPENPMX_OMEGABLOCK_MAX; i++)
		memset(&dest->omega[i], 0, sizeof(OMEGABLOCKSTYPE));

	/* Test by making a POPMODEL which may catch some errors */
	*popmodel = popmodel_init(dest->theta, dest->omega, dest->sigma, &errctx);
	if (errctx.len)
		goto failed;
	if (!args->silent) {
		fprintf(stdout, "popparam reload \"%s\"\n", filename);
		popmodel_information(0, popmodel, 0);
	}

	/* cleanup */
	/* return with success */
	if (filename_buffer)
		free(filename_buffer);
	return 0;

failed:
	/* failure path */
	free(filename_buffer);

	if (!args->silent) {
		if (errctx.len)
			fprintf(stderr, "%s", errctx.errmsg);
	}

	/* if it failed but it was optional, then dont pass error */
	if (!args->optional) 
		fatal(0, "%s: reload failed but not optional\n", __func__);

	return 1;
}

/// For the configuration object there are various settings.
///
/// - `.filename="...",` The filename is used as an .ext file to
/// load the population parameters (theta, omega, and sigma values)
/// and save them in the OPENPMX object. The default behavior is to load
/// the population paramaters from the filename of the destination
/// OPENPMX object with a .ext extension added.
/// - `.force=true,` The loaded population paramaters are copied over
/// those of the destinantion ignoring any mismatch in structure. The
/// default is to exit the program if the structures differ.
/// - `.optional=true,` If reloading fails then the destination OPENPMX
/// object is not modified and the program continues. The default
/// behavior is to consider this a fatal error and the program exits.
/// - `.silent=true,` Setting this supresses a message showing the
/// structure of the newly loaded population paramaters.

void pmx_reload_popparam(OPENPMX* dest, RELOADCONFIG* args)
{
	POPMODEL popmodel;
	let err = reload_popparam(dest, args, &popmodel);

	if (!err) {
/// A successful reload invalidates the objective function result.
		dest->result = (PMXRESULT) { 0 };

/// A successful reload invalidates any existing individual level data
/// in the destinantion OPENPMX object because the numer of etas could
/// have been changed. 
		var state = dest->state;
		if (state) {
			/* destruct and reconstruct the individual level data. We have
			 * to do it locally and then memcpy it over because IDATA has
			 * const members. */
			let idata = &state->idata;
			let advanfuncs = state->advanfuncs;
			idata_destruct(idata);
			var newidata = idata_construct(&advanfuncs->recordinfo,
										   popmodel.ntheta,
										   popmodel.nomega,
										   popmodel.nsigma,
										   advanfuncs->nstate,
										   advanfuncs->advanconfig->imodelfields.size,
										   advanfuncs->advanconfig->predictfields.size);
			memcpy(idata, &newidata, sizeof(newidata));
		}
	}
}

