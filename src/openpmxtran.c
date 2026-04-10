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

/// This file implements a translation program for structured control
/// files into compilable C source to perform model
/// estimation/simulation etc.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <unistd.h>
#include <stdarg.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <string.h>

#include "../include/openpmx.h"
#include "utils/c22.h"
#include "utils/vector.c"
#include "utils/errctx.c"
#include "utils/getdelim.c"

#define OPENPMX_FILE_RECORDS_MAX	10

static int strip_enclosing(char* s, const char terminal_begin, const char terminal_end)
{
	let len = strlen(s);
	if (len == 0)
		return 1;

	/* check at end */
	char* a = s + len - 1;
	if (*a != terminal_end)
		return 1;

	/* check at begin */
	char* b = s;
	if (*b != terminal_begin)
		return 1;

	*a = 0;
	++b;

	memmove(s, b, strlen(b) + 1);
	return 0;
}

static void strip_comments_robust(char *s)
{
    typedef enum { 
        NORMAL, 
        IN_STRING, 
        IN_CHARACTER,
        IN_COMMENT_SINGLE, 
        IN_COMMENT_MULTI,
        POTENTIAL_COMMENT_START,
        POTENTIAL_COMMENT_END
    } State;

    State state = NORMAL;
    char *src = s;
    char *dst = s;

    while (*src) {
        switch (state) {
            case NORMAL:
                if (*src == '"') {
                    state = IN_STRING;
                    *dst++ = *src;
                } else if (*src == '\'') {
                    state = IN_CHARACTER;
                    *dst++ = *src;
                } else if (*src == '/') {
                    state = POTENTIAL_COMMENT_START;
                } else {
                    *dst++ = *src;
                }
                break;

            case IN_STRING:
                if (*src == '\\' && *(src + 1) != '\0') {
                    // Handle any escaped character in string: \" or \\ etc.
                    *dst++ = *src++;
                    *dst++ = *src;
                } else if (*src == '"') {
                    state = NORMAL;
                    *dst++ = *src;
                } else {
                    *dst++ = *src;
                }
                break;

            case IN_CHARACTER:
                if (*src == '\\' && *(src + 1) != '\0') {
                    // Handle any escaped character in literal: \' or \\ etc.
                    *dst++ = *src++;
                    *dst++ = *src;
                } else if (*src == '\'') {
                    state = NORMAL;
                    *dst++ = *src;
                } else {
                    *dst++ = *src;
                }
                break;

            case POTENTIAL_COMMENT_START:
                if (*src == '/') {
                    state = IN_COMMENT_SINGLE;
                } else if (*src == '*') {
                    state = IN_COMMENT_MULTI;
                } else {
                    // False alarm: restore the slash and process current char
                    state = NORMAL;
                    *dst++ = '/';
                    // We stay on the same src char to re-evaluate it in NORMAL
                    continue; 
                }
                break;

            case IN_COMMENT_SINGLE:
                if (*src == '\n') {
                    state = NORMAL;
                    *dst++ = '\n'; // Keep newline for line numbering
                }
                break;

            case IN_COMMENT_MULTI:
                if (*src == '*') {
                    state = POTENTIAL_COMMENT_END;
                } else if (*src == '\n') {
                    *dst++ = '\n'; // Keep newline for line numbering
                }
                break;

            case POTENTIAL_COMMENT_END:
                if (*src == '/') {
                    state = NORMAL;
                } else if (*src == '*') {
                    // Sequence like ***/ - stay in POTENTIAL_END
                } else {
                    state = IN_COMMENT_MULTI;
                    if (*src == '\n') *dst++ = '\n';
                }
                break;
        }
        src++;
    }
    
    // Clean up if the file ends on a trailing slash
    if (state == POTENTIAL_COMMENT_START) *dst++ = '/';
    
    *dst = '\0';
}

static VECPTR split_sections(char *buf)
{
    VECPTR sections = { 0 };

    if (buf[0] == '$')
        vector_append(sections, 0);    // synthetic empty preamble

    vector_append(sections, buf);      // preamble text, or first $ section
    for (char *p = buf; *p != '\0'; p++) {
        if (*p == '\n' && *(p + 1) == '$') {
            *p = '\0';
            vector_append(sections, p + 1);
            p++;
        }
    }
    return sections;
}

typedef struct {
	char* args;		/* points to things in () at the start  */
	char* content;	/* points to after ')' */
} section_content;

static section_content extract_args_content(char *payload, ERRCTX* errctx)
{
    // '(' must appear as first char in payload
    char *open = payload;
    if (*open != '(') {
		errctx_add(errctx, "%s: missing '('\n", __func__);
        return (section_content) { 0 };
	}

	// Walk forward to find the closing ')', error on any nested '('
	char *close = open + 1;
	while (*close && *close != ')' && *close != '(')
		close++;
	if (*close == '(') {
		errctx_add(errctx, "%s: unexpected nested '('\n", __func__);
		return (section_content) { 0 };
	}
	if (*close == '\0') {
		errctx_add(errctx, "%s: missing ')'\n", __func__);
		return (section_content) { 0 };
	}

    *open  = '\0';      // terminate keyword
    *close = '\0';      // terminate args

	strip_firstlast_space(open + 1);
	strip_firstlast_space(close + 1);

    return (section_content) {
		.args = open + 1,
		.content = close + 1,
	};
}

typedef typeof(((OPENPMX){0}).omega[0].type) OMEGA_TYPE;

typedef struct {
	const char* filename;
	struct {
		const char* file;
		const char* code;
	} data;
	struct {
		const char* method;
		const char* init;
	} advan;
	struct {
		VECPTR args;
		const char* code;
	} imodel;
	struct {
		VECPTR args;
		const char* code;
	} predict;
	struct {
		const char* code;
	} diffeqn;
	struct {
		const char* code;
		bool nonmem_init;
	} theta;
	struct {
		OMEGA_TYPE type;
		VECPTR args;
		const char* code;
	} omega[OPENPMX_OMEGABLOCK_MAX];
	int nomega;
	struct {
		VECPTR args;
	} sigma;
	struct {
		const char* code;
	} main;
	struct {
		const char* name;
		const char* content;
	} files[OPENPMX_FILE_RECORDS_MAX];
	int nfiles;
} RESULT;

static void result_destroy(RESULT* result)
{
	vector_free(result->imodel.args);
	vector_free(result->predict.args);
	forcount(i, result->nomega)
		vector_free(result->omega[i].args);
	vector_free(result->sigma.args);
}

/// Any text preceeding the sections defined by $... is ignored and can
/// be used for a description of the control file.
static void handle_preamble(RESULT* result, char* payload, ERRCTX* errctx) 
{ 
	(void)result;
	(void)payload;
	(void)errctx; /* nothing to fail */
}

/// The `$DATA(...)`	 section defines the data filename and any code for
/// manipulation of the data file before analysis.
static void handle_data(RESULT* result, char* payload, ERRCTX* errctx) 
{ 
	if (result->data.file) {
		errctx_add(errctx, "%s: second $DATA()\n", __func__);
		return;
	}

	let text = extract_args_content(payload, errctx);
	if (errctx->len)
		return;

	strip_enclosing(text.args, '"', '"');
	result->data.file = text.args;
	result->data.code = text.content;
}

/// The `$ADVAN(...)` section defines the advancer type and any options
/// for the advancer.
static void handle_advan(RESULT* result, char* payload, ERRCTX* errctx) 
{
	if (result->advan.method) {
		errctx_add(errctx, "%s: second $ADVAN()\n", __func__);
		return;
	}
	
	let text = extract_args_content(payload, errctx);
	if (errctx->len)
		return;

	result->advan.method = text.args;
	result->advan.init = text.content;
}

/// The `$IMODEL(...)` section defines the model paramaters and the code
/// to initialize them.
static void handle_imodel(RESULT* result, char* payload, ERRCTX* errctx) 
{ 
	if (result->imodel.code) {
		errctx_add(errctx, "%s: %s: second $IMODEL()\n", __func__);
		return;
	}
	
	let text = extract_args_content(payload, errctx);
	if (errctx->len)
		return;
	
	let err = get_delim_tokens(text.args, &result->imodel.args, GET_DELIM_SEP_ANY);
	if (err) {
		errctx_add(errctx, "%s: %s: parse error\n", __func__, payload);
		return;
	}
	result->imodel.code = text.content;
}

/// The `$PREDICT(...)` section defines the prediction paramaters and the
/// code to initialize them. The predition for model fitting is returned
/// as `Y`.
static void handle_predict(RESULT* result, char* payload, ERRCTX* errctx) 
{
	if (result->predict.code) {
		errctx_add(errctx, "%s: %s: second $PREDICT()\n", __func__);
		return;
	}
	
	let text = extract_args_content(payload, errctx);
	if (errctx->len)
		return;

	let err = get_delim_tokens(text.args, &result->predict.args, GET_DELIM_SEP_ANY);
	if (err) {
		errctx_add(errctx, "%s: %s: parse error\n", __func__, payload);
		return;
	}
	result->predict.code = text.content;
}

/// The `$DIFFEQN` section defines the differential equation for the
/// advan solvers that require this.
static void handle_diffeqn(RESULT* result, char* payload, ERRCTX* errctx) 
{ 
	if (result->diffeqn.code) {
		errctx_add(errctx, "%s: %s: second $DIFFEQN()\n", __func__);
		return;
	}

	strip_firstlast_space(payload);
	result->diffeqn.code = payload;
}

/// The `$THETA` section defines the fixed effects. These can be either
/// given as `{ lower, initial, upper, ESTIMATE },` or in NONMEM style
/// as `(lower, initial, upper)` or as `(value FIXED)`.
static void handle_theta(RESULT* result, char* payload, ERRCTX* errctx) 
{ 
	if (result->theta.code) {
		errctx_add(errctx, "%s: %s: second $THETA\n", __func__);
		return;
	}

	strip_firstlast_space(payload);

	result->theta.nonmem_init = false;
	if (strchr(payload, '{') == 0) 
		result->theta.nonmem_init = true;
	
	result->theta.code = payload;
}

static void _handle_omega_helper(RESULT* result,
								 char* payload,
								 const OMEGA_TYPE omega_type,
								 ERRCTX* errctx) 
{
	let n = result->nomega;
	if (n >= OPENPMX_OMEGABLOCK_MAX) {
		errctx_add(errctx, "%s: too many omega blocks\n", __func__);
		return;
	}

	/* allow the direct C initialization of an omega block */
	if (payload[0] != '(') {
		result->omega[n].type = OMEGA_INVALID;
		vector_resize(result->omega[n].args, 0);
		result->omega[n].code = payload;
		goto block_done;
	}

	/* the other "(...)" style */
	let text = extract_args_content(payload, errctx);
	if (errctx->len)
		return;

	result->omega[n].type = omega_type;
	let err = get_delim_tokens(text.args, &result->omega[n].args, GET_DELIM_SEP_ANY);
	if (err) {
		errctx_add(errctx, "%s: parse error \"%s\"\n", __func__, payload);
		return;
	}
	if (strlen(text.content)) {
		errctx_add(errctx, "%s: parse error \"%s\"\n", __func__, text.content);
		return;
	}
	result->omega[n].code = 0;

block_done:
	result->nomega = n + 1;
}

/// The `$OMEGA()` section defines a diagonal omega matrix block of random
/// effects. The values are separated by `,` or whitespace. Negative
/// variances are treated as fixed. 
static void handle_omega(RESULT* result, char* payload, ERRCTX* errctx) 
{
	/* if there is a direct C-code initializer, then that C-code is
	 * output for this block and OMEGA_DIAG is overridden */
	
	_handle_omega_helper(result, payload, OMEGA_DIAG, errctx); 
}

/// The `$OMEGABLOCK()` section defines a block omega matrix block of
/// random effects. The values are separated by `,` or whitespace.
/// Negative variances along the diagonal are treated as fixed. 
static void handle_omegablock(RESULT* result, char* payload, ERRCTX* errctx) 
{
	if (payload[0] != '(') {
		errctx_add(errctx, "%s: $OMEGABLOCK not followed by '('\n", __func__);
		return; 
	}
	
	_handle_omega_helper(result, payload, OMEGA_BLOCK, errctx); 
}

/// The `$OMEGASAME(...)` section defines a block which copies the values
/// from lower indicies in the omega matrix. This is used for
/// inter-occasion variability, i.e. multiple etas with variances and
/// covariances matching those of lower indicies.
static void handle_omegasame(RESULT* result, char* payload, ERRCTX* errctx) 
{ 
	if (payload[0] != '(') {
		errctx_add(errctx, "%s: $OMEGASAME not followed by '('\n", __func__);
		return; 
	}
	
	_handle_omega_helper(result, payload, OMEGA_SAME, errctx); 
}

/// The `$SIGMA()` section defines a diagonal omega matrix block of
/// residual variability. The values are separated by `,` or whitespace.
/// Negative variances are treated as fixed. 
static void handle_sigma(RESULT* result, char* payload, ERRCTX* errctx) 
{
	char* input;

	/* allow the direct C initialization of an omega block */
	if (payload[0] != '(') {
		input = payload;

	/* the other "(...)" style */
	} else {
		let text = extract_args_content(payload, errctx);
		if (errctx->len)
			return;
		input = text.args;

		if (strlen(text.content)) {
			errctx_add(errctx, "%s: parse error \"%s\"\n", __func__, text.content);
			return;
		}
	}

	/* get sigma array */
	let err = get_delim_tokens(input, &result->sigma.args, GET_DELIM_SEP_ANY);
	if (err) {
		errctx_add(errctx, "%s: parse error \"%s\"\n", __func__, payload);
		return;
	}
}
/// The `$MAIN` section defines C code which will run after problem
/// initialization. 
static void handle_main(RESULT* result, char* payload, ERRCTX* errctx) 
{ 
	if (result->main.code) {
		errctx_add(errctx, "%s: second $MAIN\n", __func__);
		return;
	}

	strip_firstlast_space(payload);
	result->main.code = payload;
}

/// Sections defined by `$FILE(...)` are written to the control filename
/// appended by '.' and then the name given. For the `graphics.R` script
/// this can be used to provide a configuration file for modifying the
/// graphics produced.
static void handle_file(RESULT* result, char* payload, ERRCTX* errctx) 
{ 
	let text = extract_args_content(payload, errctx);
	if (errctx->len)
		return;

	strip_enclosing(text.args, '"', '"');
	let n = result->nfiles;
	result->files[n].name = text.args;
	result->files[n].content = text.content;
	++result->nfiles;
}

// ---------------------------------------------------------------------------
// Extract the keyword from a section string.
// A section looks like "$THETA\n   ..." or "$DATA(\"file.csv\")\n..."
// We copy characters up to the first whitespace, '(', or '\0'.
// ---------------------------------------------------------------------------
static void section_keyword(const char *section, char *out, int out_size)
{
    int i = 0;
    while (i < out_size - 1 && section[i] != '\0'
                            && section[i] != '('
                            && section[i] != ' '
                            && section[i] != '\t'
                            && section[i] != '\n'
                            && section[i] != '\r') {
        out[i] = section[i];
        i++;
    }
    out[i] = '\0';
}

/* fill in RESULT from sections */
static RESULT result_construct(VECPTR* sections, ERRCTX* errctx)
{
	/* dispatch table */
	static struct {
		const char name[16];          
		void (*handler)(RESULT* result, char* payload, ERRCTX* ctx);
	} dispatch_table[] = {
		{ "",             handle_preamble },   // section[0]: text before first $
		{ "$DATA",        handle_data     },
		{ "$ADVAN",       handle_advan    },
		{ "$IMODEL",      handle_imodel   },
		{ "$PREDICT",     handle_predict  },
		{ "$DIFFEQN",     handle_diffeqn  },
		{ "$THETA",       handle_theta    },
		{ "$OMEGA",       handle_omega    },
		{ "$OMEGABLOCK",  handle_omegablock },
		{ "$OMEGASAME",   handle_omegasame },
		{ "$SIGMA",       handle_sigma    },
		{ "$MAIN",        handle_main     },
		{ "$FILE",        handle_file     },
	};

	RESULT ret = { 0 };
	forvector(i, *sections) {
		let section = sections->mutptr[i];

		// Section 0 is the preamble — no $ keyword, use the NULL-name handler
		if (i == 0) {
			dispatch_table[0].handler(&ret, section, errctx);
			continue;
		}

		char keyword[sizeof(dispatch_table[0].name)];
		section_keyword(section, keyword, sizeof(keyword));

		// Walk the table (skip index 0, that is the preamble catch-all)
		int matched = 0;
		for (int h = 1; h < (int)countof(dispatch_table); h++) {
			let name = dispatch_table[h].name;
			if (streq(keyword, name)) {
				let payload = section + strlen(keyword);
				dispatch_table[h].handler(&ret, payload, errctx);
				if (errctx->len) 
					goto failed;

				matched = 1;
				break;
			}
		}
		if (!matched) {
			errctx_add(errctx, "%s: unknown section '%s'\n", __func__, keyword);
			goto failed;
		}
	}
	return ret;

failed:
	return (RESULT) { 0 };
}

typedef struct {
	char* mutptr;
	int size;
} file_content;

static void file_content_destroy(file_content* file)
{
	free(file->mutptr);
}

static file_content read_file(const char* filename, ERRCTX* errctx)
{
	char* buffer = 0;
	var fp = fopen(filename, "rb");
	if (!fp) {
		errctx_add(errctx, "%s: fopen \"%s\" failed\n", __func__, filename);
		goto failed;
	}

    if (fseek(fp, 0, SEEK_END) != 0) {
		errctx_add(errctx, "%s: fseek end failed\n", __func__);
		goto failed;
	}

    long filesize = ftell(fp);
	if (filesize == -1L) {
		errctx_add(errctx, "%s: ftell failed\n", __func__);
		goto failed;
	}
    
    if (fseek(fp, 0, SEEK_SET) != 0) {
		errctx_add(errctx, "%s: fseek begin failed\n", __func__);
		goto failed;
	}
	
	buffer = mallocvar(char, filesize + 1);
    if (!buffer) {
		errctx_add(errctx, "%s: malloc failed\n", __func__);
		goto failed;
	}

    let bytesRead = fread(buffer, 1, filesize, fp);
    if (ferror(fp)) {
		errctx_add(errctx, "%s: fread failed\n", __func__);
		goto failed;
	}
    buffer[bytesRead] = '\0';

    if (bytesRead == 0) {
		errctx_add(errctx, "%s: file empty\n", __func__);
		goto failed;
	}

    if (buffer[bytesRead - 1] != '\n') {
		errctx_add(errctx, "%s: file must end with newline\n", __func__);
		goto failed;
	}

	fclose(fp);
	return (file_content) {
		.mutptr = buffer,
		.size = bytesRead,
	};

failed:
	if (buffer)
		free(buffer);
	if (fp)
		fclose(fp);
	return (file_content) { 0 };
}

typedef struct {
	VECPTR header;
	VECPTR elems;
} DATAINFO;

static void datainfo_destroy(DATAINFO* datainfo)
{
	vector_free(datainfo->header);
	vector_free(datainfo->elems);
}

static double checked_atof(const char* str, ERRCTX* errctx)
{
	char *end;
	errno = 0;
	let val = strtod(str, &end);
	if (errno != 0 || end == str || *end != '\0') {
		errctx_add(errctx, "%s: invalid data \"%s\"\n", __func__, str);
		return 0;
	}
	return val;
}

static DATAINFO datainfo_create(file_content* dataptr, ERRCTX* errctx)
{
	DATAINFO ret = { 0 };
	char* line = dataptr->mutptr;
	VECPTR v = { 0 };

	/* sanity check */
	if (!line || !*line) {
		errctx_add(errctx, "%s: data file is empty\n", __func__);
		return ret;
	}

	/* get first line, terminate it so we can read tokens */
	char* next_newline = strchr(line, '\n');
	if (next_newline) 
		*next_newline = '\0';

	/* determine separator based on header line */
	var sep = GET_DELIM_SEP_WHITESPACE;
	if (strchr(line, ','))
		sep = GET_DELIM_SEP_COMMA;

	/* parse the header */
	strip_firstlast_space(line);
	let err = get_delim_tokens(line, &ret.header, sep);
	if (err || ret.header.size == 0) {
		errctx_add(errctx, "%s: parse error or empty data header\n", __func__);
		goto failed;
	}

	/* data row loop */
	int linenum = 1;
	line = next_newline ? next_newline + 1 : NULL;

	vector_reserve(v, ret.header.size);
	while (line && *line) {
		linenum++;
		next_newline = strchr(line, '\n');
		if (next_newline) 
			*next_newline = '\0';

		vector_resize(v, 0);
		strip_firstlast_space(line);
		get_delim_tokens(line, &v, sep);

		/* strict size check - fails on blank lines or count mismatch */
		if (v.size != ret.header.size) {
			errctx_add(errctx, "%s: data row %i has %i elements, header has %i\n",
					   __func__, linenum, v.size, ret.header.size);
			goto failed;
		}

		/* check for valid values */
		forvector_val(value, v) {
			if (!streq(value, ".")) {
				(void)checked_atof(value, errctx);
				if (errctx->len) 
					goto failed;
			}
		}
		
		vector_appendn(ret.elems, v.mutptr, v.size);

		/* advance to start of next line or stop if buffer ends */
		line = next_newline ? next_newline + 1 : NULL;
	}

	/* final check - ensure at least one row of data exists */
	if (linenum < 2) {
		errctx_add(errctx, "%s: no data rows after header\n", __func__);
		goto failed;
	}
	vector_free(v);
	return ret;

failed:
	vector_free(v);
	vector_free(ret.header);
	vector_free(ret.elems);
	return (DATAINFO) { 0 };
}

static int checked_atoi(const char* str, ERRCTX* errctx)
{
	char *end;
	errno = 0;
	let val = strtol(str, &end, 10);
	if (errno != 0 || end == str || *end != '\0') {
		errctx_add(errctx, "%s: invalid value \"%s\"\n", __func__, str);
		return 0;
	}
	return (int)val;
}

static void write_code(FILE* fp, const char* code, const char* missing)
{
	if (code && strlen(code))
		fprintf(fp, "\t%s\n", code);
	else 
		fprintf(fp, "%s\n", missing);
}

static void write_args(FILE* fp, const char* pre, const VECPTR* args, const char* post)
{
	forvector_val(v, *args) 
		fprintf(fp, "\t%s%s = %s%s;\n", pre, v, post, v);
}

static void template_filename(FILE* fp, 
							const RESULT* const result, 
							const DATAINFO* const datainfo, 
							ERRCTX* errctx)
{
	(void)datainfo;
	(void)errctx;
	fprintf(fp, "\"%s\"\n", result->filename);
}

static void template_data_array(FILE* fp, 
								const RESULT* const result, 
								const DATAINFO* const datainfo, 
								ERRCTX* errctx)
{
	(void)result;
	(void)errctx;
	
	fprintf(fp, "typedef struct RECORD {\n");
	let header = &datainfo->header;
	forvector_val(v, *header)
		fprintf(fp, "\tdouble %s;\n", v);
	fprintf(fp, "} RECORD;\n");
	fprintf(fp, "RECORD OPENPMXTRAN_DATA_NAME[] = {\n");

	int ncols = 0;
	let elems = &datainfo->elems;
	forvector_val(v, *elems) {
		if (ncols == 0)
			fprintf(fp, "\t{ ");

/// Within the datafile a `.` as a data entry is replaced by NAN.
		if (streq(v, ".")) {
			v = "NAN";

		} else {
			(void)checked_atof(v, errctx);
			if (errctx->len)
				return;
		}

		fprintf(fp, "%s, ", v);
		
		++ncols;
		if (ncols == header->size) {
			fprintf(fp, "},\n");
			ncols = 0;
		}
	}
	fprintf(fp, "};\n");
}

static void template_imodel_fields_declare(FILE* fp, 
										const RESULT* const result, 
										const DATAINFO* const datainfo, 
										ERRCTX* errctx)
{
	(void)datainfo;
	(void)errctx;
	forvector_val(v, result->imodel.args) 
		fprintf(fp, "\tdouble %s;\n", v);
}

static void template_imodel_fields_declare_set(FILE* fp, 
											   const RESULT* const result, 
											   const DATAINFO* const datainfo, 
											   ERRCTX* errctx)
{
	(void)datainfo;
	(void)errctx;
	write_args(fp, "double ", &result->imodel.args, "_imodel->");
}

static void template_imodel_code(FILE* fp, 
								 const RESULT* const result, 
								 const DATAINFO* const datainfo, 
								 ERRCTX* errctx)
{
	(void)datainfo;

	/* its an error to have no model code */
	if (!result->imodel.code || strlen(result->imodel.code) == 0) {
		errctx_add(errctx, "%s: no $IMODEL code\n", __func__);
		return;
	}

	write_code(fp, result->imodel.code, "\t/* no imodel code */\n");

	/* suppress warnings that no theta was used if there are no thetas
	 * available */
	if (!result->theta.code || strlen(result->theta.code) == 0)
		fprintf(fp, "\t(void)_theta; /* no thetas */\n");

	/* suppress warnings that no omega was used if there are no omegas
	 * available */
	if (result->nomega == 0)
		fprintf(fp, "\t(void)_eta; /* no omegas */\n");
}

static void template_imodel_fields_set(FILE* fp, 
									   const RESULT* const result, 
									   const DATAINFO* const datainfo, 
									   ERRCTX* errctx)
{
	(void)datainfo;
	(void)errctx;
	
	write_args(fp, "_imodel->", &result->imodel.args, "");
}

static void template_diffeqn_code(FILE* fp, 
								  const RESULT* const result, 
								  const DATAINFO* const datainfo, 
								  ERRCTX* errctx)
{
	(void)errctx;
	(void)datainfo;
	
	write_code(fp, result->diffeqn.code, 
			   "\t/* no diffeqn code */\n\t(void)_dadt;\n");
}

static void template_record_fields_define(FILE* fp, 
										  const RESULT* const result, 
										  const DATAINFO* const datainfo, 
										  ERRCTX* errctx)
{
	(void)result;
	(void)errctx;
	
	forvector_val(v, datainfo->header) 
		fprintf(fp, "\tconst double %s = _record->%s; (void)%s;\n", v, v, v);
}

static void template_imodel_fields_define(FILE* fp, 
										  const RESULT* const result, 
										  const DATAINFO* const datainfo, 
										  ERRCTX* errctx)
{
	(void)datainfo;
	(void)errctx;
	
	forvector_val(v, result->imodel.args) 
		fprintf(fp, "\tconst double %s = _imodel->%s; (void)%s;\n", v, v, v);
}

static void template_predparams_fields_declare(FILE* fp, 
											const RESULT* const result, 
											const DATAINFO* const datainfo, 
											ERRCTX* errctx)
{
	(void)datainfo;
	(void)errctx;
	
	forvector_val(v, result->predict.args) 
		fprintf(fp, "\tdouble %s;\n", v);
}

static void template_predict_code(FILE* fp, 
								  const RESULT* const result, 
								  const DATAINFO* const datainfo, 
								  ERRCTX* errctx)
{
	(void)errctx;
	(void)datainfo;
	
	write_code(fp, result->predict.code, 
			   "\t/* no predict code */\n\t(void)_err;\n");
}

static void template_predparams_fields_set(FILE* fp, 
										   const RESULT* const result, 
										   const DATAINFO* const datainfo, 
										   ERRCTX* errctx)
{
	(void)datainfo;
	(void)errctx;
	
	write_args(fp, "_predparams->", &result->predict.args, "");
}

static void template_record_fields_writeable_define(FILE* fp, 
													const RESULT* const result, 
													const DATAINFO* const datainfo, 
													ERRCTX* errctx)
{
	(void)result;
	(void)errctx;
	
	write_args(fp, "double ", &datainfo->header, "_record->");
}

static void template_data_preprocess_code(FILE* fp, 
										  const RESULT* const result, 
										  const DATAINFO* const datainfo, 
										  ERRCTX* errctx)
{
	(void)errctx;
	(void)datainfo;
	
	write_code(fp, result->data.code, "\t/* no data preprocess code */\n");
}

static void template_record_fields_writeback(FILE* fp, 
											 const RESULT* const result, 
											 const DATAINFO* const datainfo, 
											 ERRCTX* errctx)
{
	(void)result;
	(void)errctx;
	
	write_args(fp, "_record->", &datainfo->header, "");
}

static void template_data_config(FILE* fp, 
								 const RESULT* const result, 
								 const DATAINFO* const datainfo, 
								 ERRCTX* errctx)
{
	(void)result;
	(void)errctx;
	
	fprintf(fp, "\t\t.writeable = OPENPMXTRAN_DATA_NAME,\n");
	fprintf(fp, "\t\t.records = OPENPMXTRAN_DATA_NAME,\n");
	fprintf(fp, "\t\t.nrecords = sizeof(OPENPMXTRAN_DATA_NAME)/sizeof(OPENPMXTRAN_DATA_NAME[0]),\n");
	fprintf(fp, "\t\t._offset1 = true,\n");
	fprintf(fp, "\t\t.recordfields = {\n");
	fprintf(fp, "\t\t\t.size = sizeof(RECORD),\n");
	fprintf(fp, "\t\t\t.field = {\n");
	forvector_val(v, datainfo->header) {
		fprintf(fp, "\t\t\t\t{ .name=\"%s\", .offset=offsetof(RECORD, %s) },\n", v, v);
	}
	fprintf(fp, "\t\t\t},\n");
	fprintf(fp, "\t\t},\n");
}

static void write_fieldinfo(FILE* fp, 
							const VECPTR* args, 
							const char* type)
{
	forvector_val(v, *args) 
		fprintf(fp, "\t\t\t\t{ .name=\"%s\", .offset = offsetof(%s, %s) },\n", v, type, v);
}

static void template_imodel_fieldinfo(FILE* fp, 
									  const RESULT* const result, 
									  const DATAINFO* const datainfo, 
									  ERRCTX* errctx)
{
	(void)errctx;
	(void)datainfo;
	write_fieldinfo(fp, &result->imodel.args, "IMODEL");
}

static void template_predparams_fieldinfo(FILE* fp, 
										  const RESULT* const result, 
										  const DATAINFO* const datainfo, 
										  ERRCTX* errctx)
{
	(void)errctx;
	(void)datainfo;
	write_fieldinfo(fp, &result->predict.args, "PREDICTVARS");
}

static void template_advan_init(FILE* fp, 
								const RESULT* const result, 
								const DATAINFO* const datainfo, 
								ERRCTX* errctx)
{
	(void)datainfo;

/// Within the `$ADVAN(...)` section an advancer must be indicated.

	static const struct {
		const char* method;
		const char* funcname;
	} advan_table[] = {
		/// + `pred` calls `pmx_advan_pred()` a simple predictor
		{ "pred",                 "pmx_advan_pred"                 },
		/// + `onecomp` calls `pmx_advan_onecomp()` a one-compartment model.
		{ "onecomp",              "pmx_advan_onecomp"              },
		/// + `onecomp_depot` calls `pmx_advan_onecomp_depot()` a one-compartment model with a depot compartment.
		{ "onecomp_depot",        "pmx_advan_onecomp_depot"        },
		/// + `twocomp` calls `pmx_advan_twocomp()` a two-compartment mammilary model.
		{ "twocomp",              "pmx_advan_twocomp"              },
		/// + `threecomp` calls `pmx_advan_threecomp()` a three-compartment mammilary model.
		{ "threecomp",            "pmx_advan_threecomp"            },
		{ "diffeqn_test",         "pmx_advan_diffeqn_test"         },
		/// + `diffeqn_libgsl` calls `pmx_advan_diffeqn_libgsl()` a ODE solver from LibGSL.
		{ "diffeqn_libgsl",       "pmx_advan_diffeqn_libgsl"       },
		/// + `eigen` calls `pmx_advan_eigen()` a linear eigensystem solver. In the $IMODEL() function
		/// the eigensystem matrix must be specified by SYSMAT().
		{ "eigen",                "pmx_advan_eigen"                },
		/// + `eigen_threecomp` calls `pmx_advan_eigen_threecomp()` a linear eigensystem solver specialized
		/// to a three compartment model. The eigensystem matrix does not have to be set, it is set automatically.
		{ "eigen_threecomp",      "pmx_advan_eigen_threecomp"      },
		/// + `eigen_twocomp` calls `pmx_advan_eigen_twocomp()` a linear eigensystem solver specialized
		/// to a two compartment model. The eigensystem matrix does not have to be set, it is set automatically.
		{ "eigen_twocomp",        "pmx_advan_eigen_twocomp"        },
		/// + `eigen_onecomp_absorb` calls `pmx_advan_eigen_onecomp_absorb()` a linear eigensystem solver specialized
		/// to a one compartment model with absorbtion. The eigensystem matrix does not have to be set, it is set automatically.
		{ "eigen_onecomp_absorb", "pmx_advan_eigen_onecomp_absorb" },
	};

	let method = result->advan.method;
	const char* methodname = NULL;
	for (int i = 0; i < (int)countof(advan_table); i++) {
		if (streq(method, advan_table[i].method)) {
			methodname = advan_table[i].funcname;
			break;
		}
	}
	if (!methodname) {
		errctx_add(errctx, "%s: invalid advan \"%s\"\n", __func__, method);
		return;
	}

	fprintf(fp, "\t\t.method = %s,\n", methodname);
	fprintf(fp, "\t\t%s\n", result->advan.init);
}

static const char *my_strcasestr(const char *haystack, const char *needle)
{
	size_t needle_len = strlen(needle);
	if (needle_len == 0)
		return haystack;
	for (; *haystack; haystack++) {
		const char *h = haystack;
		const char *n = needle;
		while (*n && tolower((unsigned char)*h) == tolower((unsigned char)*n)) {
			h++;
			n++;
		}
		if (*n == '\0')
			return haystack;
	}
	return NULL;
}

static bool keyword_present(const char *s, const char *kw)
{
    const char *hit = my_strcasestr(s, kw);
    if (!hit)
		return false;
    char after = hit[strlen(kw)];
    return !isalnum((unsigned char)after) && after != '_';
}

/**
 * Parses NONMEM style (val FIXED) or (val, val, val). Returns false on
 * syntax errors like consecutive commas. Outputs to the provided FILE*.
 * The input line is modified in place (null terminators are inserted to
 *  delimit numeric tokens).
 */
static void transform_pmx_line(FILE *fp, char *line, ERRCTX* errctx)
{
	char *p = line;

	while (*p != '\0') {
		// 1. Skip leading whitespace to find the start of a group
		while (isspace((unsigned char)*p))
			p++;
		if (*p == '\0')
			break;

		// 2. Strict Check: If it's not an opening paren, it's an illegal character/comment
		if (*p != '(') {
			errctx_add(errctx, "%s: Unexpected character '%c' found outside parentheses\n", __func__, *p);
			return;
		}

		// --- Parsing inside the (...) 
		char *cursor = p + 1;
		const char *val_starts[3] = {NULL};
		int val_lens[3] = {0};
		int found = 0;
		bool is_fixed = (keyword_present(cursor, "FIXED") || keyword_present(cursor, "FIX"));

		for (int i = 0; i < 3; i++) {
			while (isspace((unsigned char)*cursor))
				cursor++;
			
			// Check for illegal characters/comments inside the parens
			if (*cursor != '-' && *cursor != '.' && !isdigit((unsigned char)*cursor) && 
				*cursor != ',' && *cursor != ')' && !isalpha((unsigned char)*cursor)) {
				errctx_add(errctx, "%s: Illegal character '%c' inside group\n", __func__, *cursor);
				return;
			}

			if (*cursor == ',') {
				errctx_add(errctx, "%s: consecutive or leading comma detected\n", __func__);
				return;
			}

			char *token_start = cursor;
			char *endptr;
			errno = 0;
			(void)strtod(cursor, &endptr);
			if (cursor == endptr) {
				errctx_add(errctx, "%s: no digit parsed at \"%s\"\n", __func__, cursor);
				break;
			}
			if (errno == ERANGE) {
				errctx_add(errctx, "%s: out of double range \"%s\"\n", __func__, cursor);
				break;
			}

			val_starts[found] = token_start;
			val_lens[found] = (int)(endptr - token_start);
			found++;

			cursor = endptr;
			while (isspace((unsigned char)*cursor))
				cursor++;
			if (*cursor == ',')
				cursor++;
			else if (*cursor == ')' || isalpha((unsigned char)*cursor))
				break;
		}
		if (found == 0) {
			errctx_add(errctx, "%s: no values found\n", __func__);
			return;
		}

		// Output logic...
		if (found > 0) {
			fprintf(fp, "    {%.*s, %.*s, %.*s, %s },\n",
					val_lens[0], val_starts[0],
					(found >= 2) ? val_lens[1] : val_lens[0], (found >= 2) ? val_starts[1] : val_starts[0],
					(found >= 3) ? val_lens[2] : val_lens[0], (found >= 3) ? val_starts[2] : val_starts[0],
					is_fixed ? "FIXED" : "ESTIMATE");
		}

		// 3. Close the group and update p
		char *close = strchr(cursor, ')');
		if (!close) {
			errctx_add(errctx, "%s: Missing closing parenthesis\n", __func__);
			return;
		}
		p = close + 1; // Move p to just after the ')' for the next loop iteration
	}
}

static void template_theta_init(FILE* fp, 
								const RESULT* const result, 
								const DATAINFO* const datainfo, 
								ERRCTX* errctx)
{
	(void)datainfo;

	/* NONMEM style initializer */
	if (result->theta.nonmem_init) {
		var codecopy = strdup(result->theta.code);
		transform_pmx_line(fp, codecopy, errctx);
		free(codecopy);

		if (errctx->len)
			return;

	/* C-style initializer */
	} else {
		write_code(fp, result->theta.code, "\t/* no theta init code */\n");
	}
}

static void _omega_args_to_code(FILE* fp, const int type, const VECPTR* args, ERRCTX* errctx)
{
	var ndim = -1;
	var print_values = true;
	var type_name = "";

	if (type == OMEGA_DIAG) {
		type_name = "OMEGA_DIAG";
		ndim = args->size;
		
	} else if (type == OMEGA_BLOCK) {
		type_name = "OMEGA_BLOCK";
		/* calculate block dimensions */
		/* https://www.wolframalpha.com/input?i2d=true&i=solve+k%5C%2844%29+n%3Dk*Divide%5B%5C%2840%29k%2B1%5C%2841%29%2C2%5D */
		let n = args->size;
		ndim = (int)floor((sqrt(8 * n + 1) - 1) / 2);

		/* double check that the size is consistant */
		let nvals = ndim * (ndim + 1) / 2;
		if (n != nvals) {
			errctx_add(errctx, "%s: omega block size inconsistant (%i,%i,%i)\n",
					   __func__, n, ndim, nvals);
			return;
		}
		
	} else if (type == OMEGA_SAME) {
		type_name = "OMEGA_SAME";
		if (args->size != 1) {
			errctx_add(errctx, "%s: omega same size is not 1 (%i)\n", __func__, args->size);
			return;
		}
		
		/* extract omega size catching errors */
		ndim = checked_atoi(args->ptr[0], errctx);
		if (errctx->len) {
			errctx_add(errctx, "%s: invalid omega same \"%s\"\n", __func__, args->ptr[0]);
			return;
		}
		print_values = false;
		
	} else {
		errctx_add(errctx, "%s: invalid omega type (%i)\n", __func__, type);
		return;
	}

	fprintf(fp, "\t\t{ %s, %i, { ", type_name, ndim);
	if (print_values) {
		forvector(j, *args) {
			let v = args->ptr[j];
			fprintf(fp, "%s, ", v);
		}
	}
	fprintf(fp, " } },\n");
}

static void template_omega_init(FILE* fp, 
								const RESULT* const result, 
								const DATAINFO* const datainfo, 
								ERRCTX* errctx)
{
	(void)datainfo;
	
	forcount(i, result->nomega) {
		var type = result->omega[i].type;
		let args = &result->omega[i].args;
		let code = result->omega[i].code;
		if (code && strlen(code)) {
			fprintf(fp, "%s", code);

		} else {
			_omega_args_to_code(fp, type, args, errctx);
			if (errctx->len)
				return;
		}
	}
}

static void template_sigma_init(FILE* fp, 
								const RESULT* const result, 
								const DATAINFO* const datainfo, 
								ERRCTX* errctx)
{
	(void)errctx;
	(void)datainfo;

	fprintf(fp, "\t\t");
	forvector_val(v, result->sigma.args)
		fprintf(fp, "%s, ", v);
	fprintf(fp, "\n");
}

static void template_main_code(FILE* fp, 
							   const RESULT* const result, 
							   const DATAINFO* const datainfo, 
							   ERRCTX* errctx)
{
	(void)errctx;
	(void)datainfo;
	write_code(fp, result->main.code, "\t/* no main code */\n");
}

#if defined(__MSYS__) || defined(__MINGW32__)
#ifndef NAME_MAX
#define NAME_MAX 255
#endif
#endif /* defined(__MSYS__) || defined(__MINGW32__) */

static void expand_template(const char* grfilename, 
							char* template, 
							const RESULT* const result, 
							const DATAINFO* const datainfo, 
							ERRCTX* errctx)
{
	static struct {
		const char name[36];
		void (*handler)(FILE* fp, 
						const RESULT* const result, 
						const DATAINFO* const datainfo, 
						ERRCTX* errctx);
	} dispatch_table[] = {
		{ "{{DEFAULT_FILENAME}}",				template_filename },
		{ "{{DATA_ARRAY}}",						template_data_array },
		{ "{{IMODEL_FIELDS_DECLARE}}",			template_imodel_fields_declare },
		{ "{{IMODEL_FIELDS_DECLARE_AND_SET}}",	template_imodel_fields_declare_set },
		{ "{{IMODEL_CODE}}",					template_imodel_code },
		{ "{{IMODEL_FIELDS_SET}}",				template_imodel_fields_set },
		{ "{{IMODEL_DIFFEQN_CODE}}",			template_diffeqn_code },
		{ "{{RECORD_FIELDS_DEFINE}}",			template_record_fields_define },
		{ "{{IMODEL_FIELDS_DEFINE}}",			template_imodel_fields_define },
		{ "{{PREDPARAMS_FIELDS_DECLARE}}", 		template_predparams_fields_declare },
		{ "{{PREDICT_CODE}}", 					template_predict_code },
		{ "{{PREDPARAMS_FIELDS_SET}}",			template_predparams_fields_set },
		{ "{{RECORD_FIELDS_WRITABLE_DEFINE}}",	template_record_fields_writeable_define },
		{ "{{DATA_PREPROCESS_CODE}}", 			template_data_preprocess_code },
		{ "{{RECORD_FIELDS_WRITEBACK}}", 		template_record_fields_writeback },
		{ "{{DATA_CONFIG}}", 					template_data_config },
		{ "{{IMODEL_FIELDINFO}}", 				template_imodel_fieldinfo },
		{ "{{PREDPARAMS_FIELDINFO}}", 			template_predparams_fieldinfo },
		{ "{{ADVAN_INIT}}", 					template_advan_init },
		{ "{{THETA_INIT}}", 					template_theta_init },
		{ "{{OMEGA_INIT}}", 					template_omega_init },
		{ "{{SIGMA_INIT}}", 					template_sigma_init },
		{ "{{MAIN_CODE}}", 						template_main_code },
	};

	/* figure out the filename */
	char c_filename[NAME_MAX + PATH_MAX];
	snprintf(c_filename, sizeof(c_filename), "%s.c", grfilename);
	var fp = fopen(c_filename, "w");
	if (!fp) {
		errctx_add(errctx, "%s: could not open \"%s\"\n", __func__, c_filename);
		goto done;
	}
	
	/* process template line by line */
	char *saveptr;
	let delim = "\n";
	char *line = strtok_r(template, delim, &saveptr); 
	while (line != NULL) {
		
		/* match a template line */
		char keyword[sizeof(dispatch_table[0].name)];
		section_keyword(line, keyword, sizeof(keyword));

		/* look for a handler for the template line, call it */
		int matched = 0;
		for (int h = 0; h < (int)countof(dispatch_table); h++) {
			let name = dispatch_table[h].name;
			if (streq(keyword, name)) {
				dispatch_table[h].handler(fp, result, datainfo, errctx);
				if (errctx->len)
					goto done;
				matched = 1;
				break;
			}
		}
		/* no match, so we just output the template */
		if (!matched) 
			fprintf(fp, "%s\n", line);

		/* go to next template line */
		line = strtok_r(NULL, delim, &saveptr);
	}

done:
	/* cleanup */
	if (fp)
		fclose(fp);

	/* if error cleanup the output C file */
	if (errctx->len)
		unlink(c_filename);
}

/* validate variable names */
int is_valid_identifier(const char *s)
{
    if (s == NULL || *s == '\0')
        return 0;

    // First character: letter or underscore
    if (!(isalpha((unsigned char)*s) || *s == '_'))
        return 0;

    s++;

    // Remaining characters: letter, digit, or underscore
    while (*s) {
        if (!(isalnum((unsigned char)*s) || *s == '_'))
            return 0;
        s++;
    }

    return 1;
}

static const char *keywords[] =
{
    "auto","break","case","char","const","continue","default","do","double",
    "else","enum","extern","float","for","goto","if","inline","int","long",
    "register","restrict","return","short","signed","sizeof","static","struct",
    "switch","typedef","union","unsigned","void","volatile","while","_Alignas",
    "_Alignof","_Atomic","_Bool","_Complex","_Generic","_Imaginary","_Noreturn",
    "_Static_assert","_Thread_local"
};

static int is_keyword(const char *s)
{
    size_t n = sizeof(keywords)/sizeof(keywords[0]);
    for (size_t i = 0; i < n; i++) {
        if (strcmp(s, keywords[i]) == 0)
            return 1;
    }
    return 0;
}

static int cmp_str(const void *a, const void *b)
{
	const char *sa = *(const char **)a;
	const char *sb = *(const char **)b;
	return strcmp(sa, sb);
}

typedef VECTOR(const char*) CONSTVECPTR;

static void check_names(const RESULT* const result, const DATAINFO* const datainfo, ERRCTX* errctx)
{
	CONSTVECPTR allnames = { 0 };

	/* put all names into a list and sort them */
	vector_appendn(allnames, datainfo->header.ptr, datainfo->header.size);
	vector_appendn(allnames, result->imodel.args.ptr, result->imodel.args.size);
	vector_appendn(allnames, result->predict.args.ptr, result->predict.args.size);
    qsort(allnames.ptr, allnames.size, sizeof(allnames.ptr[0]), cmp_str);

	/* doubles will be sequential, easy to detect */
    for (int i = 1; i < allnames.size; i++) {
		let v1 = allnames.ptr[i - 1];
		let v2 = allnames.ptr[i];
		if (strcmp(v1, v2) == 0) {
			errctx_add(errctx, "%s: name \"%s\" has multiple definitions\n", __func__, v1);
			goto done;
		}
	}

	/* look for invalid names */
	forvector_val(v, allnames) {
		if (!is_valid_identifier(v)) {
			errctx_add(errctx, "%s: invalid name \"%s\"\n", __func__, v);
			goto done;
		}
		if (is_keyword(v)) {
			errctx_add(errctx, "%s: \"%s\" is a C keyword\n", __func__, v);
			goto done;
		}
	}

	/* look for reserved names everywhere */
    forvector_val(v, allnames) {
		if ((strcmp(v, "Y") == 0) ||
			(strcmp(v, "YHAT") == 0) ||
			(strcmp(v, "YHATVAR") == 0) ||
			(strcmp(v, "PRED") == 0) ||
			(strcmp(v, "OBJ") == 0) ||
			(strcmp(v, "INEVAL") == 0) ||
			(strcmp(v, "NAN") == 0)) {
			errctx_add(errctx, "%s: name \"%s\" is reserved\n", __func__, v);
			goto done;
		}
	}
done:
    vector_free(allnames);
}

static void write_files(const RESULT* result, ERRCTX* errctx)
{
#define FILENAME (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
	
	char fname[PATH_MAX];
	forcount(i, result->nfiles) {
		snprintf(fname, sizeof(fname), "%s.%s", result->filename, result->files[i].name);
		printf("%s: write \"%s\"\n", FILENAME, fname);
		
		var f = fopen(fname, "w");
		if (!f) {
			errctx_add(errctx, "%s: unable to open file \"%s\"\n", __func__, fname);
			continue;
		}
		
		fputs(result->files[i].content, f);
		fputc('\n', f);
		if (ferror(f)) 
			errctx_add(errctx, "%s: error writing to file \"%s\"\n", __func__, fname);

		fclose(f);
	}
}

int main(int argc, char* argv[])
{
	if (argc != 2 && argc != 3) {
		fprintf(stderr, "Usage: %s <controlfile> [datafile]\n", argv[0]);
		exit(EXIT_FAILURE);
	}
 
	let filename = strrchr("/" __FILE__, '/') + 1;
	let grfilename = argv[1];
	var datafilename = (const char*)0;
	if (argc == 3 && strlen(argv[2]))
		datafilename = argv[2];

	/* read in control file */
	ERRCTX errctx = { };
	var controltext = read_file(grfilename, &errctx);
	if (errctx.len) {
		fprintf(stderr, "%s: %s", filename, errctx.errmsg);
		exit(EXIT_FAILURE);
	}
	strip_comments_robust(controltext.mutptr);

	/* split into sections with $... */
	var sections = split_sections(controltext.mutptr);
	
	/* handle each section, write into result */
	var result = result_construct(&sections, &errctx);
	if (errctx.len) {
		fprintf(stderr, "%s: %s", filename, errctx.errmsg);
		exit(EXIT_FAILURE);
	}
	result.filename = grfilename;

	/* make sure data file is not given twice */
	if (result.data.file && datafilename) {
		fprintf(stderr, "%s: data file in $DATA() and also on command line\n", filename);
		exit(EXIT_FAILURE);
	}
	var datafile_to_use = datafilename;
	if (!datafile_to_use)
		datafile_to_use = result.data.file;
	if (!datafile_to_use) {
		fprintf(stderr, "%s: no data file\n", filename);
		exit(EXIT_FAILURE);
	}

	/* read in datafile and parse it */
	var datatext = read_file(datafile_to_use, &errctx);
	if (errctx.len) {
		fprintf(stderr, "%s: %s", filename, errctx.errmsg);
		exit(EXIT_FAILURE);
	}
	var datainfo = datainfo_create(&datatext, &errctx);
	if (errctx.len) {
		fprintf(stderr, "%s: %s", filename, errctx.errmsg);
		exit(EXIT_FAILURE);
	}

	/* do some error checking before we go further */
	check_names(&result, &datainfo, &errctx);
	if (errctx.len) {
		fprintf(stderr, "%s: %s", filename, errctx.errmsg);
		exit(EXIT_FAILURE);
	}
	
	/* process the template */
	extern char openpmxtran_template[];
	expand_template(grfilename, openpmxtran_template, &result, &datainfo, &errctx);
	if (errctx.len) {
		fprintf(stderr, "%s: %s", filename, errctx.errmsg);
		exit(EXIT_FAILURE);
	}

	/* write the files */
	write_files(&result, &errctx);
	if (errctx.len) {
		fprintf(stderr, "%s: %s", filename, errctx.errmsg);
		exit(EXIT_FAILURE);
	}

	/* cleanup */
	datainfo_destroy(&datainfo);
	result_destroy(&result);
	file_content_destroy(&datatext);
	file_content_destroy(&controltext);
	vector_free(sections);
	exit(EXIT_SUCCESS);
}

char openpmxtran_template[] =
"/*\n"
" * Code below this message is generated by OpenPMX (https://github.com/deleveld/openpmx).\n"
" * Copyright (c) 2024 Douglas Eleveld.\n"
" *\n"
" * This program is free software: you can redistribute it and/or modify\n"
" * it under the terms of the GNU General Public License as published by\n"
" * the Free Software Foundation, version 3.\n"
" *\n"
" * This program is distributed in the hope that it will be useful, but\n"
" * WITHOUT ANY WARRANTY; without even the implied warranty of\n"
" * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU\n"
" * General Public License for more details.\n"
" *\n"
" * You should have received a copy of the GNU General Public License\n"
" * along with this program. If not, see <http://www.gnu.org/licenses/>.\n"
" */\n"
"\n"
"#include <math.h>\n"
"#include <stdlib.h>\n"
"#include <stddef.h>\n"
"#include <assert.h>\n"
"#include <string.h>\n"
"//#define _GNU_SOURCE\n"
"//#include <fenv.h>\n"
"#include \"openpmx.h\"\n"
"\n"
"#define OPENPMXTRAN_DEFAULT_FILENAME \\\n"
"{{DEFAULT_FILENAME}}\n"
"\n"
"{{DATA_ARRAY}}\n"
"\n"
"typedef struct IMODEL {\n"
"{{IMODEL_FIELDS_DECLARE}}\n"
"} IMODEL;\n"
"\n"
"static void imodel_init(IMODEL* const _imodel,\n"
"						 ADVANSTATE* const _advanstate)\n"
"{\n"
"	const double* _theta = _advanstate->current.popparam->theta;\n"
"	const double* _eta = _advanstate->current.popparam->eta;\n"
"	const RECORD* const _record = _advanstate->current.record;\n"
"\n"
"	/* allow access to RECORD fields */\n"
"{{RECORD_FIELDS_DEFINE}}\n"
"\n"
"	/* declare and set IMODEL fields */\n"
"{{IMODEL_FIELDS_DECLARE_AND_SET}}\n"
"\n"
"#define THETA(i) 		((const double)_theta[(i)-1])\n"
"#define ETA(i) 		((const double)_eta[(i)-1])\n"
"#define A(i) 			((const double)_advanstate->current.state[(i)-1])\n"
"#define ALAG(i,t) 		pmx_advan_amtlag(_advanstate, (i)-1, (t))\n"
"#define BIOAVAIL(i,f)	pmx_advan_bioavail(_advanstate, (i)-1, (f))\n"
"#define STATETIME 		((const double)_advanstate->current.statetime)\n"
"#define INITTIME(t) 	pmx_advan_inittime(_advanstate, (t))\n"
"#define A_0(i,v) 		pmx_advan_state_init(_advanstate, (i)-1, (v))\n"
"#define SYSMAT(...)	pmx_advan_eigen_sysmat(_advanstate, (double[]){ __VA_ARGS__ })\n"
"\n"
"{{IMODEL_CODE}}\n"
"\n"
"#undef THETA\n"
"#undef ETA\n"
"#undef A\n"
"#undef ALAG\n"
"#undef BIOAVAIL\n"
"#undef STATETIME\n"
"#undef INITTIME\n"
"#undef A_0\n"
"#undef SYSMAT\n"
"\n"
"	/* set IMODEL fields */\n"
"{{IMODEL_FIELDS_SET}}\n"
"}\n"
"\n"
"typedef struct PREDICTVARS {\n"
"{{PREDPARAMS_FIELDS_DECLARE}}\n"
"} PREDICTVARS;\n"
"\n"
"static void imodel_diffeqn(double _dadt[],\n"
"						   const IMODEL* const _imodel,\n"
"						   const RECORD* const _record,\n"
"						   const double* const _state,\n"
"						   const POPPARAM* const _popparam,\n"
"						   const double T)\n"
"{\n"
"	const double* const _theta = _popparam->theta;\n"
"	const double* const _eta = _popparam->eta;\n"
"	(void)_record;\n"
"	(void)_theta;\n"
"	(void)_eta;\n"
"	(void)_state;\n"
"	(void)T;\n"
"\n"
"	/* allow access to RECORD fields */\n"
"{{RECORD_FIELDS_DEFINE}}\n"
"\n"
"	/* declare IMODEL fields */\n"
"{{IMODEL_FIELDS_DEFINE}}\n"
"\n"
"#define THETA(i) 	((const double)_theta[(i)-1])\n"
"#define ETA(i) 	((const double)_eta[(i)-1])\n"
"#define A(i) 		((const double)_state[(i)-1])\n"
"#define DADT(i) 	(_dadt[(i)-1])\n"
"\n"
"{{IMODEL_DIFFEQN_CODE}}\n"
"\n"
"#undef THETA\n"
"#undef ETA\n"
"#undef A\n"
"#undef DADT\n"
"}\n"
"\n"
"static double imodel_predict(const IMODEL* const _imodel,\n"
"							 const PREDICTSTATE* const _current,\n"
"							 const double* const _err,\n"
"							 PREDICTVARS* _predparams)\n"
"{\n"
"	const RECORD* const _record = _current->record;\n"
"	const POPPARAM* const _popparam = _current->popparam;\n"
"	(void) _record;\n"
"	(void) _predparams;\n"
"	(void) _popparam;\n"
"	double Y = NAN;\n"
"\n"
"	/* allow access to RECORD fields */\n"
"{{RECORD_FIELDS_DEFINE}}\n"
"\n"
"	/* allow access to IMODEL fields */\n"
"{{IMODEL_FIELDS_DEFINE}}\n"
"\n"
"	/* declare PREDICTVARS fields */\n"
"{{PREDPARAMS_FIELDS_DECLARE}}\n"
"\n"
"#define THETA(i) 	((const double)_popparam->theta[(i)-1])\n"
"#define ETA(i) 	((const double)_popparam->eta[(i)-1])\n"
"#define ERR(i) 	((const double)_err[(i)-1])\n"
"#define EPS(i) 	((const double)_err[(i)-1])\n"
"#define A(i) 		((const double)_current->state[(i)-1])\n"
"\n"
"{{PREDICT_CODE}}\n"
"\n"
"#undef THETA\n"
"#undef ETA\n"
"#undef A\n"
"#undef ERR\n"
"#undef EPS\n"
"\n"
"	/* set PREDICTVARS fields */\n"
"{{PREDPARAMS_FIELDS_SET}}\n"
"\n"
"	return Y;\n"
"}\n"
"\n"
"static void openpmxtran_data_preprocess_callback(RECORD* _record)\n"
"{\n"
"	/* allow write access to RECORD fields */\n"
"{{RECORD_FIELDS_WRITABLE_DEFINE}}\n"
"\n"
"	/* allow user to change data */\n"
"\tdo {\n"
"{{DATA_PREPROCESS_CODE}}\n"
"\t} while (0);\n"
"\n"
"	/* write back RECORD fields */\n"
"{{RECORD_FIELDS_WRITEBACK}}\n"
"\n"
"	return;\n"
"}\n"
"\n"
"static OPENPMX openpmx = (OPENPMX) {\n"
"	.filename = OPENPMXTRAN_DEFAULT_FILENAME,\n"
"	.nthread = 0,\n"
"	.data = { \n"
"{{DATA_CONFIG}}\n"
"	},\n"
"	.advan = {\n"
"		.init = imodel_init,\n"
"		.predict = imodel_predict,\n"
"		.diffeqn = imodel_diffeqn,\n"
"		.imodelfields = {\n"
"			.size = sizeof(IMODEL),\n"
"			.field = {\n"
"{{IMODEL_FIELDINFO}}\n"
"			},\n"
"		},\n"
"		.predictfields = {\n"
"			.size = sizeof(PREDICTVARS),\n"
"			.field = {\n"
"{{PREDPARAMS_FIELDINFO}}\n"
"			},\n"
"		},\n"
"{{ADVAN_INIT}}\n"
"	},\n"
"	.theta = {\n"
"{{THETA_INIT}}\n"
"	},\n"
"	.omega = {\n"
"{{OMEGA_INIT}}\n"
"	},\n"
"	.sigma = {\n"
"{{SIGMA_INIT}}\n"
"	},\n"
"};\n"
"\n"
"void estimate_config(ESTIMCONFIG* estimconfig)\n"
"{\n"
"	pmx_estimate(&openpmx, estimconfig);\n"
"}\n"
"#define estimate(...) pmx_estimate(&openpmx, &(ESTIMCONFIG){ __VA_ARGS__ })\n"
"\n"
"void fastestimate_config(ESTIMCONFIG* estimconfig)\n"
"{\n"
"	pmx_fastestimate(&openpmx, estimconfig);\n"
"}\n"
"#define fastestimate(...) pmx_fastestimate(&openpmx, &(ESTIMCONFIG){ __VA_ARGS__ })\n"
"\n"
"void evaluate_config(STAGE1CONFIG* stage1config)\n"
"{\n"
"	pmx_evaluate(&openpmx, stage1config);\n"
"}\n"
"#define evaluate(...) pmx_evaluate(&openpmx, &(STAGE1CONFIG){ __VA_ARGS__ })\n"
"\n"
"void simulate_config(const SIMCONFIG* simconfig)\n"
"{\n"
"	pmx_simulate(&openpmx, simconfig);\n"
"}\n"
"#define simulate(...) pmx_simulate(&openpmx, &(SIMCONFIG){ __VA_ARGS__ })\n"
"\n"
"void predict(void)\n"
"{\n"
"	pmx_predict(&openpmx);\n"
"}\n"
"\n"
"void predict_pred(void)\n"
"{\n"
"	pmx_predict_pred(&openpmx);\n"
"}\n"
"\n"
"#define table(f, ...)	pmx_table(&openpmx, (f), &(TABLECONFIG){ __VA_ARGS__ })\n"
"\n"
"#define reload(...) pmx_reload_popparam(&openpmx, &(RELOADCONFIG){ __VA_ARGS__ })\n"
"\n"
"#define set_theta(i, ...) pmx_set_theta(&openpmx, (i), &(typeof(((OPENPMX){0}).theta[0])){ __VA_ARGS__ })\n"
"\n"
"static void openpmxtran_data_preprocess(void)\n"
"{\n"
"	RECORD* data = openpmx.data.writeable;\n"
"	if (!data)\n"
"		return;\n"
"	const int nrecords = openpmx.data.nrecords;\n"
"	for (int i=0; i<nrecords; i++) {\n"
"		RECORD r = data[i];\n"
"		openpmxtran_data_preprocess_callback(&r);\n"
"		data[i] = r;\n"
"	}\n"
"}\n"
"\n"
"extern void server_queue(void);\n"
"\n"
"int main(void)\n"
"{\n"
"	/* enabling exceptions can sometimes help debugging */\n"
"	/* feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW); */\n"
"\n"
"	/* do any preprocessing defined in the $DATA code */\n"
"	openpmxtran_data_preprocess();\n"
"\n"
"	/* wait in the server queue if necessary */\n"
"	server_queue();\n"
"\n"
"{{MAIN_CODE}}\n"
"\n"
"	pmx_cleanup(&openpmx);\n"
"	return EXIT_SUCCESS;\n"
"}\n";

