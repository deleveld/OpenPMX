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

#include "../include/openpmx.h"
#include "utils/c22.h"
#include "utils/vector.c"
#include "utils/errctx.c"

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
	memmove(s, a, strlen(a) + 1);
}

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
                } else if (*src == '/') {
                    state = POTENTIAL_COMMENT_START;
                } else {
                    *dst++ = *src;
                }
                break;

            case IN_STRING:
                if (*src == '\\' && *(src + 1) == '"') {
                    // Handle escaped quotes: \"
                    *dst++ = *src++;
                    *dst++ = *src;
                } else if (*src == '"') {
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
                    state = NORMAL;
                    *dst++ = '/'; // It wasn't a comment, restore the slash
                    *dst++ = *src;
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
                    // Stay in this state (handles sequences like ***/)
                } else {
                    state = IN_COMMENT_MULTI;
                    if (*src == '\n') *dst++ = '\n';
                }
                break;
        }
        src++;
    }
    
    // Finalize: handle edge cases where file ends mid-potential-comment
    if (state == POTENTIAL_COMMENT_START) *dst++ = '/';
    *dst = '\0';
}

typedef VECTOR(char*) VECPTR;

typedef enum {
	GET_DELIM_SEP_COMMA,		/* force comma separator */
	GET_DELIM_SEP_WHITESPACE,	/* force whitespace separator */
	GET_DELIM_SEP_ANY,			/* use comma separator if comma exists */
} GET_DELIM_SEP;

static int get_delim_tokens_comma(char* line, VECPTR* namevec)
{
	if (line[0] == ',')
		return 1;
	
	char *running = line;
	char *token;
	while ((token = strsep(&running, ",")) != NULL) {
		strip_firstlast_space(token);
		if (strlen(token)) {
			vector_append(*namevec, token);

		/* missing token */
		} else
			return 1;
	}
	return 0;
}

static void get_delim_tokens_whitespace(char* line, VECPTR* namevec)
{
	char *saveptr;
	let delim = " \t\r\n";
	char *token = strtok_r(line, delim, &saveptr); 
	while (token != NULL) {
		vector_append(*namevec, token);
		token = strtok_r(NULL, delim, &saveptr);
	}
}

static int get_delim_tokens(char* line, VECPTR* namevec, const GET_DELIM_SEP sep)
{
	var s = sep;
	if (sep == GET_DELIM_SEP_ANY) {
		if (strchr(line, ',')) 
			s = GET_DELIM_SEP_COMMA;
		else
			s = GET_DELIM_SEP_WHITESPACE;
	}
	
	/* comma separated, means a single comma separates */
	if (s == GET_DELIM_SEP_COMMA) {
		return get_delim_tokens_comma(line, namevec);

	/* whitespace separated, means multiple whitespace separates */
	} else if (s == GET_DELIM_SEP_WHITESPACE) {
		get_delim_tokens_whitespace(line, namevec);
		
	/* this should never happen */
	} else {
		assert(0);
	}
		
	return 0;
}

static void split_sections(char *buf, VECPTR* sections)
{
    /* If the file starts with $, synthesise an empty preamble section
     * so that sections[0] is always the pre-$ preamble text. */
    if (buf[0] == '$')
        vector_append(*sections, "");

    /* sections[0] always begins at buf (may be the synthesised "" above
     * if buf[0]=='$', or real preamble text otherwise). */
    vector_append(*sections, buf);

    /* Walk the buffer looking for a newline immediately followed by '$'.
     * When found, terminate the current section at the newline and record
     * a pointer to the '$' as the start of the next section. */
    for (char *p = buf; *p != '\0'; p++) {
        if (*p == '\n' && *(p + 1) == '$') {
            *p = '\0';                        /* terminate current section */
            vector_append(*sections, p + 1); /* next section starts at '$' */
            p++;                              /* skip the '\0'; loop p++ will
                                                then step onto the char after
                                                '$', which is correct */
        }
    }
}

typedef struct {
	char* args;
	char* content;
} section_content;

static section_content extract_section_args(char *payload, ERRCTX* errctx)
{
    // '(' must appear as first char in payload
    char *open = payload;
    if (*open != '(') {
		add_errctx(errctx, "%s: missing '('\n", __func__);
        return (section_content) { 0 };
	}

    // Walk forward tracking nesting depth to find the matching ')'
    int depth = 0;
    char *p = open;
    char *close = NULL;
    while (*p) {
        if      (*p == '(') depth++;
        else if (*p == ')') {
            depth--;
            if (depth == 0) { close = p; break; }
        }
        p++;
    }
    if (!close) {
		add_errctx(errctx, "%s: missing ')'\n", __func__);
        return (section_content) { 0 };
	}

    *open  = '\0';      // terminate keyword
    *close = '\0';      // terminate args

    return (section_content) {
		.args = open + 1,
		.content = close + 1,
	};
}

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
	} theta;
	struct {
		const char* type;
		VECPTR args;
	} omega[OPENPMX_OMEGABLOCK_MAX];
	int nomega;
	struct {
		VECPTR args;
	} sigma;
	struct {
		const char* code;
	} main;
} RESULT;

static void handle_preamble(RESULT* result, const int keywordlen, char* section, ERRCTX* errctx) 
{ 
	(void)result;
	(void)keywordlen;
	(void)section;
	(void)errctx; /* nothing to fail */
}

static void handle_data(RESULT* result, const int keywordlen, char* section, ERRCTX* errctx) 
{ 
	let payload = section + keywordlen;
	let text = extract_section_args(payload, errctx);
	if (errctx->len)
		return;

	strip_firstlast_space(text.args);
	strip_enclosing(text.args, '"', '"');
	result->data.file = text.args;

	strip_firstlast_space(text.content);
	result->data.code = text.content;
}

static void handle_advan(RESULT* result, const int keywordlen, char* section, ERRCTX* errctx) 
{ 
	let payload = section + keywordlen;
	let text = extract_section_args(payload, errctx);
	if (errctx->len)
		return;

	strip_firstlast_space(text.args);
	strip_firstlast_space(text.content);

	result->advan.method = text.args;
	result->advan.init = text.content;
}

static void handle_imodel(RESULT* result, const int keywordlen, char* section, ERRCTX* errctx) 
{ 
	let payload = section + keywordlen;
	let text = extract_section_args(payload, errctx);
	if (errctx->len)
		return;
	
	strip_firstlast_space(text.args);
	strip_firstlast_space(text.content);

	let err = get_delim_tokens(text.args, &result->imodel.args, GET_DELIM_SEP_ANY);
	if (err) {
		add_errctx(errctx, "%s: %s: parse error\n", __func__, section);
		return;
	}
	result->imodel.code = text.content;
}

static void handle_predict(RESULT* result, const int keywordlen, char* section, ERRCTX* errctx) 
{ 
	let payload = section + keywordlen;
	let text = extract_section_args(payload, errctx);
	if (errctx->len)
		return;

	strip_firstlast_space(text.args);
	strip_firstlast_space(text.content);

	let err = get_delim_tokens(text.args, &result->predict.args, GET_DELIM_SEP_ANY);
	if (err) {
		add_errctx(errctx, "%s: %s: parse error\n", __func__, section);
		return;
	}
	result->predict.code = text.content;
}

static void handle_diffeqn(RESULT* result, const int keywordlen, char* section, ERRCTX* errctx) 
{ 
	(void)errctx; /* nothing to fail */

	let payload = section + keywordlen;
	strip_firstlast_space(payload);
	result->diffeqn.code = payload;
}

static void handle_theta(RESULT* result, const int keywordlen, char* section, ERRCTX* errctx) 
{ 
	(void)errctx; /* nothing to fail */

	let payload = section + keywordlen;
	strip_firstlast_space(payload);
	result->theta.code = payload;
}

static void handle_omega(RESULT* result, const int keywordlen, char* section, ERRCTX* errctx) 
{ 
	let payload = section + keywordlen;
	let text = extract_section_args(payload, errctx);
	if (errctx->len)
		return;

	strip_firstlast_space(text.args);
	strip_firstlast_space(text.content);

	let n = result->nomega;
	if (n >= OPENPMX_OMEGABLOCK_MAX) {
		add_errctx(errctx, "%s: too many omega blocks\n", __func__);
		return;
	}

	result->omega[n].type = section; /* truncated by extract_section_args()... as intended */
	let err = get_delim_tokens(text.args, &result->omega[n].args, GET_DELIM_SEP_ANY);
	if (err) {
		add_errctx(errctx, "%s: parse error \"%s\"\n", __func__, section);
		return;
	}
	result->nomega = n + 1;

	if (strlen(text.content)) {
		add_errctx(errctx, "%s: parse error \"%s\"\n", __func__, text.content);
		return;
	}
}

static void handle_omegablock(RESULT* result, const int keywordlen, char* section, ERRCTX* errctx) 
{ 
	handle_omega(result, keywordlen, section, errctx);
}

static void handle_omegasame(RESULT* result, const int keywordlen, char* section, ERRCTX* errctx) 
{ 
	handle_omega(result, keywordlen, section, errctx);
}

static void handle_sigma(RESULT* result, const int keywordlen, char* section, ERRCTX* errctx) 
{ 
	let payload = section + keywordlen;
	let text = extract_section_args(payload, errctx);
	if (errctx->len)
		return;

	strip_firstlast_space(text.args);
	strip_firstlast_space(text.content);

	let err = get_delim_tokens(text.args, &result->sigma.args, GET_DELIM_SEP_ANY);
	if (err) {
		add_errctx(errctx, "%s: parse error \"%s\"\n", __func__, section);
		return;
	}
		
	if (strlen(text.content)) {
		add_errctx(errctx, "%s: parse error \"%s\"\n", __func__, text.content);
		return;
	}
}

static void handle_main(RESULT* result, const int keywordlen, char* section, ERRCTX* errctx) 
{ 
	(void)errctx; /* nothing to fail */

	let payload = section + keywordlen;
	strip_firstlast_space(payload);
	result->main.code = payload;
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

// ---------------------------------------------------------------------------
// Dispatch table  (NULL name = catch-all for the leading preamble section)
// ---------------------------------------------------------------------------
struct {
	const char  *name;          
	void (*handler)(RESULT* result, 
					const int keywordlen, 
					char* content, 
					ERRCTX* ctx);
} handlers[] = {
	{ NULL,           handle_preamble },   // section[0]: text before first $
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
};

// ---------------------------------------------------------------------------
// dispatch_sections()
//   sections  – pointer array from split_sections()
// ---------------------------------------------------------------------------
static void dispatch_sections(RESULT* result, VECPTR* sections, ERRCTX* errctx)
{
    forvector(i, *sections) {
		var section = sections->mutptr[i];

        // Section 0 is the preamble — no $ keyword, use the NULL-name handler
        if (i == 0) {
            handlers[0].handler(result, 0, section, errctx);
            continue;
        }

        char keyword[512];
        section_keyword(section, keyword, sizeof(keyword));

        // Walk the table (skip index 0, that is the preamble catch-all)
        int matched = 0;
        for (int h = 1; h < (int)countof(handlers); h++) {
			let name = handlers[h].name;
            if (streq(keyword, name)) {
				let keywordlen = strlen(keyword);
                handlers[h].handler(result, keywordlen, section, errctx);
                if (errctx->len) 
					return;

                matched = 1;
                break;
            }
        }
        if (!matched) {
			add_errctx(errctx, "%s: unknown section '%s'\n", __func__, keyword);
			return;
        }
    }
}

static char* read_file(const char* filename, ERRCTX* errctx)
{
	var fp = fopen(filename, "rb");
	if (!fp) {
		add_errctx(errctx, "%s: fopen \"%s\" failed\n", __func__, filename);
		return 0;
	}

    if (fseek(fp, 0, SEEK_END) != 0) {
		add_errctx(errctx, "%s: fseek end failed\n", __func__);
		fclose(fp);
		return 0;
	}
    let filesize = ftell(fp);
	if (filesize == -1L) {
		add_errctx(errctx, "%s: ftell failed\n", __func__);
		fclose(fp);
		return 0;
	}
    
    if (fseek(fp, 0, SEEK_SET) != 0) {
		add_errctx(errctx, "%s: fseek begin failed\n", __func__);
		fclose(fp);
		return 0;
	}
	
    var buffer = mallocvar(char, filesize + 1);
    if (!buffer) {
		add_errctx(errctx, "%s: malloc failed\n", __func__);
		fclose(fp);
		return 0;
	}

    let bytesRead = fread(buffer, 1, filesize, fp);
    if (ferror(fp)) {
		add_errctx(errctx, "%s: fread failed\n", __func__);
		free(buffer);
		fclose(fp);
		return 0;
	}
    buffer[bytesRead] = '\0';

	fclose(fp);
	return buffer;
}

typedef struct {
	VECPTR header;
	VECPTR elems;
} DATAINFO;

static void parse_datafile(char* datatext, DATAINFO* datainfo, ERRCTX* errctx)
{
	/* read in first line as header */
	char *saveptr;
	let delim = "\n";
	char *line = strtok_r(datatext, delim, &saveptr); 
	if (!line) {
		add_errctx(errctx, "%s: parse header failed\n", __func__);
		return;
	}
	
	/* use the same separator for the whole file */
	var sep = GET_DELIM_SEP_WHITESPACE;
	if (strchr(line, ','))
		sep = GET_DELIM_SEP_COMMA;
	
	/* parse header */
	var header = &datainfo->header;
	let err = get_delim_tokens(line, header, sep);
	if (err) {
		add_errctx(errctx, "%s: parse header tokens failed\n", __func__);
		return;
	}
	var linenum = 1;

	/* parse the rest of the lines */
	line = strtok_r(NULL, delim, &saveptr);
	while (line) {
		++linenum;

		VECPTR v = { };
		vector_reserve(v, header->size);
		let err = get_delim_tokens(line, &v, sep);
		if (err) {
			add_errctx(errctx, "%s: parse data elements failed line %i\n", __func__, linenum);
			vector_free(v);
			return;
		}
		if (v.size != header->size) {
			add_errctx(errctx, "%s: data row %i has %i elements, header has %i\n", __func__, linenum, v.size, header->size);
			vector_free(v);
			return;
		}
		
		/* accept the line, copy over the elements */
		vector_appendn(datainfo->elems, v.mutptr, v.size);
		vector_free(v);
		
		/* get next line */
		line = strtok_r(NULL, delim, &saveptr);
	}
	
	if (linenum == 1) 
		add_errctx(errctx, "%s: no data rows\n", __func__);
}

static int checked_atoi(const char* str, ERRCTX* errctx)
{
	char *end;
	errno = 0;
	let val = strtol(str, &end, 10);
	if (errno != 0 || end == str || *end != '\0') {
		add_errctx(errctx, "%s: invalid value \"%s\"\n", __func__, str);
		return 0;
	}
	return (int)val;
}

static void write_code(FILE* fp, const char* code, const char* name)
{
	if (code && strlen(code))
		fprintf(fp, "\t%s\n", code);
	else 
		fprintf(fp, "\t/* no %s code */\n", name);
}

static void write_args(FILE* fp, const char* pre, const VECPTR* args, const char* post)
{
	forvector(i, *args) {
		let v = args->ptr[i];
		fprintf(fp, "\t%s%s = %s%s;\n", pre, v, post, v);
	}
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
	(void)datainfo;
	(void)errctx;
	
	fprintf(fp, "typedef struct RECORD {\n");
	let header = &datainfo->header;
	forvector(i, *header) {
		let v = header->ptr[i];
		fprintf(fp, "\tdouble %s;\n", v);
	}
	fprintf(fp, "} RECORD;\n");
	fprintf(fp, "RECORD OPENPMXTRAN_DATA_NAME[] = {\n");

	int ncols = 0;
	let elems = &datainfo->elems;
	forvector(i, *elems) {
		if (ncols == 0)
			fprintf(fp, "\t{ ");
		
		var v = elems->ptr[i];
		if (streq(v, "."))
			v = "NAN";
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
	forvector(i, result->imodel.args) {
		let v = result->imodel.args.ptr[i];
		fprintf(fp, "\tdouble %s;\n", v);
	}
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
	(void)errctx;
	(void)datainfo;
	write_code(fp, result->imodel.code, "imodel");
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
	write_code(fp, result->diffeqn.code, "diffeqn");

	/* avoid a useless warning, but keep it if needed */
	if (result->diffeqn.code == 0)
		fprintf(fp, "\t(void)_dadt;\n");
}

static void template_record_fields_define(FILE* fp, 
										const RESULT* const result, 
										const DATAINFO* const datainfo, 
										ERRCTX* errctx)
{
	(void)result;
	(void)errctx;
	forvector(i, datainfo->header) {
		let v = datainfo->header.ptr[i];
		fprintf(fp, "\tconst double %s = _record->%s; (void)%s;\n", v, v, v);
	}
}

static void template_imodel_fields_define(FILE* fp, 
										const RESULT* const result, 
										const DATAINFO* const datainfo, 
										ERRCTX* errctx)
{
	(void)datainfo;
	(void)errctx;
	forvector(i, result->imodel.args) {
		let v = result->imodel.args.ptr[i];
		fprintf(fp, "\tconst double %s = _imodel->%s; (void)%s;\n", v, v, v);
	}
}

static void template_predparams_fields_declare(FILE* fp, 
											const RESULT* const result, 
											const DATAINFO* const datainfo, 
											ERRCTX* errctx)
{
	(void)datainfo;
	(void)errctx;
	forvector(i, result->predict.args) {
		let v = result->predict.args.ptr[i];
		fprintf(fp, "\tdouble %s;\n", v);
	}
}

static void template_predict_code(FILE* fp, 
								const RESULT* const result, 
								const DATAINFO* const datainfo, 
								ERRCTX* errctx)
{
	(void)errctx;
	(void)datainfo;
	write_code(fp, result->predict.code, "predict");
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
	write_code(fp, result->data.code, "data preprocess");
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
	forvector(i, datainfo->header) {
		let v = datainfo->header.ptr[i];
		fprintf(fp, "\t\t\t\t{ .name=\"%s\", .offset=offsetof(RECORD, %s) },\n", v, v);
	}
	fprintf(fp, "\t\t\t},\n");
	fprintf(fp, "\t\t},\n");
}

static void write_fieldinfo(FILE* fp, const VECPTR* args, const char* type)
{
	forvector(i, *args) {
		let v = args->ptr[i];
		fprintf(fp, "\t\t\t\t{ .name=\"%s\", .offset = offsetof(%s, %s) },\n", v, type, v);
	}
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
	(void)errctx;
	(void)datainfo;
	
	let method = result->advan.method;
	var methodname = "";

	/// + `pred` calls `pmx_advan_pred()` a simple predictor
	if (streq(method, "pred")) {
		methodname = "pmx_advan_pred";
		
	/// + `onecomp` calls `pmx_advan_onecomp()` a one-compartment model.
	} else if (streq(method, "onecomp")) 
		methodname = "pmx_advan_onecomp";
		
	/// + `onecomp_depot` calls `pmx_advan_onecomp_depot()` a one-compartment model with a depot compartment.
	else if (streq(method, "onecomp_depot"))
		methodname = "pmx_advan_onecomp_depot";
		
	/// + `twocomp` calls `pmx_advan_twocomp()` a two-compartment model mammilary model.
	else if (streq(method, "twocomp"))
		methodname = "pmx_advan_twocomp";
		
	/// + `threecomp` calls `pmx_advan_threecomp()` a three-compartment model mammilary model.
	else if (streq(method, "threecomp"))
		methodname = "pmx_advan_threecomp";
		
	else if (streq(method, "diffeqn_test")) {
		methodname = "pmx_advan_diffeqn_test";
		
	/// + `diffeqn` calls `pmx_advan_diffeqn_libgsl()` a ODE solver from LibGSL.
	} else if (streq(method, "diffeqn_libgsl")) {
		methodname = "pmx_advan_diffeqn_libgsl";
		
	/// + `eigen` calls `pmx_advan_eigen()` a linear eigensystem solver. In the $IMODEL() function 
	/// the eigensystem matrix must be specified by SYSMAT(). 
	} else if (streq(method, "eigen")) 
		methodname = "pmx_advan_eigen";
		
	/// + `eigen_threecomp` calls `pmx_advan_eigen_threecomp()` a linear eigensystem solver specialized
	/// to a three compartment model. The eigensystem matrix does not have to be set, it is set 
	/// automatically.
	else if (streq(method, "eigen_threecomp"))
		methodname = "pmx_advan_eigen_threecomp";
		
	/// + `eigen_twocomp` calls `pmx_advan_eigen_twocomp()` a linear eigensystem solver specialized
	/// to a two compartment model. The eigensystem matrix does not have to be set, it is set 
	/// automatically.
	else if (streq(method, "eigen_twocomp"))
		methodname = "pmx_advan_eigen_twocomp";
		
	/// + `eigen_onecomp_absorb` calls `pmx_advan_eigen_onecomp_absorb()` a linear eigensystem solver specialized
	/// to a one compartment model with absorbtion. The eigensystem matrix does not have to be set, it is set 
	/// automatically.
	else if (streq(method, "eigen_onecomp_absorb"))
		methodname = "pmx_advan_eigen_onecomp_absorb";
	
	else {
		add_errctx(errctx, "%s: invalid advan \"%s\"\n", __func__, method);
		return;
	}

	fprintf(fp, "\t\t.method = %s,\n", methodname);
	fprintf(fp, "\t\t%s\n", result->advan.init);
}

static void template_theta_init(FILE* fp, 
								const RESULT* const result, 
								const DATAINFO* const datainfo, 
								ERRCTX* errctx)
{
	(void)errctx;
	(void)datainfo;
	write_code(fp, result->theta.code, "theta init");
}

static void template_omega_init(FILE* fp, 
								const RESULT* const result, 
								const DATAINFO* const datainfo, 
								ERRCTX* errctx)
{
	(void)errctx;
	(void)datainfo;
	
	forcount(i, result->nomega) {
		var type = result->omega[i].type;
		let args = result->omega[i].args;
		var ndim = -1;

		var print_values = true;
		if (strcmp(type, "$OMEGA") == 0) {
			type = "OMEGA_DIAG";
			ndim = args.size;
			
		} else if (strcmp(type, "$OMEGABLOCK") == 0) {
			type = "OMEGA_BLOCK";
			/* calculate block dimensions */
			/* https://www.wolframalpha.com/input?i2d=true&i=solve+k%5C%2844%29+n%3Dk*Divide%5B%5C%2840%29k%2B1%5C%2841%29%2C2%5D */
			let n = args.size;
			ndim = (int)floor((sqrt(8 * n + 1) - 1) / 2);

			/* double check that the size is consistant */
			let nvals = ndim * (ndim + 1) / 2;
			if (n != nvals) {
				add_errctx(errctx, "%s: omega block size inconsistant (%i,%i,%i)\n", __func__, n, ndim, nvals);
				return;
			}
			
		} else if (strcmp(type, "$OMEGASAME") == 0) {
			type = "OMEGA_SAME";
			if (args.size != 1) {
				add_errctx(errctx, "%s: omega same size is not 1 (%i)\n", __func__, args.size);
				return;
			}
			
			/* extract omega size catching errors */
			ndim = checked_atoi(args.ptr[0], errctx);
			if (errctx->len) {
				add_errctx(errctx, "%s: invalid omega same \"%s\"\n", __func__, args.ptr[0]);
				return;
			}
			print_values = false;
			
		} else {
			add_errctx(errctx, "%s: invalid omega type \"%s\"\n", __func__, type);
			return;
		}
	
		fprintf(fp, "\t\t{ %s, %i, { ", type, ndim);
		if (print_values) {
			forvector(j, args) {
				let v = args.ptr[j];
				fprintf(fp, "%s, ", v);
			}
		}
		fprintf(fp, " } },\n");
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
	forvector(i, result->sigma.args) {
		let v = result->sigma.args.ptr[i];
		fprintf(fp, "%s, ", v);
	}
	fprintf(fp, "\n");
}

static void template_main_code(FILE* fp, 
							   const RESULT* const result, 
							   const DATAINFO* const datainfo, 
							   ERRCTX* errctx)
{
	(void)errctx;
	(void)datainfo;
	write_code(fp, result->main.code, "main");
}

struct {
	const char  *name;          
	void (*handler)(FILE* fp, 
					const RESULT* const result, 
					const DATAINFO* const datainfo, 
					ERRCTX* errctx);
} template_entries[] = {
	{ "{{OPENPMXTRAN_DEFAULT_FILENAME}}",				template_filename },
	{ "{{OPENPMXTRAN_DATA_ARRAY}}",						template_data_array },
	{ "{{OPENPMXTRAN_IMODEL_FIELDS_DECLARE}}",			template_imodel_fields_declare },
	{ "{{OPENPMXTRAN_IMODEL_FIELDS_DECLARE_AND_SET}}",	template_imodel_fields_declare_set },
	{ "{{OPENPMXTRAN_IMODEL_CODE}}",					template_imodel_code },
	{ "{{OPENPMXTRAN_IMODEL_FIELDS_SET}}",				template_imodel_fields_set },
	{ "{{OPENPMXTRAN_IMODEL_DIFFEQN_CODE}}",			template_diffeqn_code },
	{ "{{OPENPMXTRAN_RECORD_FIELDS_DEFINE}}",			template_record_fields_define },
	{ "{{OPENPMXTRAN_IMODEL_FIELDS_DEFINE}}",			template_imodel_fields_define },
	{ "{{OPENPMXTRAN_PREDPARAMS_FIELDS_DECLARE}}", 		template_predparams_fields_declare },
	{ "{{OPENPMXTRAN_PREDICT_CODE}}", 					template_predict_code },
	{ "{{OPENPMXTRAN_PREDPARAMS_FIELDS_SET}}",			template_predparams_fields_set },
	{ "{{OPENPMXTRAN_RECORD_FIELDS_WRITABLE_DEFINE}}",	template_record_fields_writeable_define },
	{ "{{OPENPMXTRAN_DATA_PREPROCESS_CODE}}", 			template_data_preprocess_code },
	{ "{{OPENPMXTRAN_RECORD_FIELDS_WRITEBACK}}", 		template_record_fields_writeback },
	{ "{{OPENPMXTRAN_DATA_CONFIG}}", 					template_data_config },
	{ "{{OPENPMXTRAN_IMODEL_FIELDINFO}}", 				template_imodel_fieldinfo },
	{ "{{OPENPMXTRAN_PREDPARAMS_FIELDINFO}}", 			template_predparams_fieldinfo },
	{ "{{OPENPMXTRAN_ADVAN_INIT}}", 					template_advan_init },
	{ "{{OPENPMXTRAN_THETA_INIT}}", 					template_theta_init },
	{ "{{OPENPMXTRAN_OMEGA_INIT}}", 					template_omega_init },
	{ "{{OPENPMXTRAN_SIGMA_INIT}}", 					template_sigma_init },
	{ "{{OPENPMXTRAN_MAIN_CODE}}", 						template_main_code },
};

static void expand_template(const char* grfilename, 
							char* template, 
							const RESULT* const result, 
							const DATAINFO* const datainfo, 
							ERRCTX* errctx)
{
	/* figure out the filename */
	let ext = ".c";
	var c_filename = malloc(strlen(grfilename) + strlen(ext) + 1);
	strcpy(c_filename, grfilename);
	strcat(c_filename, ext);
	var fp = fopen(c_filename, "w");
	if (!fp) {
		add_errctx(errctx, "%s: could not open \"%s\"\n", __func__, c_filename);
		goto done;
	}
	
	/* process template line by line */
	char *saveptr;
	let delim = "\n";
	char *line = strtok_r(template, delim, &saveptr); 
	while (line != NULL) {
		
		/* match a template line */
		char keyword[512];
		section_keyword(line, keyword, sizeof(keyword));

		/* look for a handler for the template line, call it */
		int matched = 0;
		for (int h = 0; h < (int)countof(template_entries); h++) {
			let name = template_entries[h].name;
			if (streq(keyword, name)) {
				template_entries[h].handler(fp, result, datainfo, errctx);
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
		
	free(c_filename);
}

static int cmp_str(const void *a, const void *b)
{
	const char *sa = *(const char **)a;
	const char *sb = *(const char **)b;
	return strcmp(sa, sb);
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

static int is_valid_identifier_strict(const char *s)
{
    if (!is_valid_identifier(s))
        return 0;
    return !is_keyword(s);
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
			add_errctx(errctx, "%s: name \"%s\" has multiple definitions\n", __func__, v1);
			goto done;
		}
	}

	/* look for invalid names */
	forvector(i, allnames) {
		let v = allnames.ptr[i];
		if (!is_valid_identifier(v)) {
			add_errctx(errctx, "%s: invalid name \"%s\"\n", __func__, v);
			goto done;
		}
		if (is_keyword(v)) {
			add_errctx(errctx, "%s: \"%s\" clashes with C keyword\n", __func__, v);
			goto done;
		}
	}

	/* look for reserved names everwhere */
    forvector(i, allnames) {
		let v = allnames.ptr[i];
		if ((strcmp(v, "Y") == 0) ||
			(strcmp(v, "YHAT") == 0) ||
			(strcmp(v, "YHATVAR") == 0) ||
			(strcmp(v, "PRED") == 0) ||
			(strcmp(v, "OBJ") == 0) ||
			(strcmp(v, "INEVAL") == 0) ||
			(strcmp(v, "NAN") == 0)) {
			add_errctx(errctx, "%s: name \"%s\" is reserved\n", __func__, v);
			goto done;
		}
	}
done:
    vector_free(allnames);
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
	strip_comments_robust(controltext);

	/* split into sections with $... */
	VECPTR sections = { };
	split_sections(controltext, &sections);
	
	/* handle each section, write into result */
	RESULT result = { 
		.filename = grfilename,
	};
	dispatch_sections(&result, &sections, &errctx);
	if (errctx.len) {
		fprintf(stderr, "%s: %s", filename, errctx.errmsg);
		exit(EXIT_FAILURE);
	}

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
	DATAINFO datainfo = { };
	parse_datafile(datatext, &datainfo, &errctx);
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

	/* cleanup datainfo */
	vector_free(datainfo.header);
	vector_free(datainfo.elems);

	/* cleanup result */
	vector_free(result.imodel.args);
	vector_free(result.predict.args);
	forcount(i, result.nomega)
		vector_free(result.omega[i].args);
	vector_free(result.sigma.args);

	/* cleanup */
	free(datatext);
	vector_free(sections);
	free(controltext);
	exit(EXIT_SUCCESS);
}

/// The variables and code defined in the block is placed into a
/// template of an OpenPMX analysis.
char openpmxtran_template[] =
"/*\n"
" * Code below this message is generated by openpmx (https://github.com/deleveld/openpmx).\n"
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
"{{OPENPMXTRAN_DEFAULT_FILENAME}}\n"
"\n"
"{{OPENPMXTRAN_DATA_ARRAY}}\n"
"\n"
"typedef struct IMODEL {\n"
"{{OPENPMXTRAN_IMODEL_FIELDS_DECLARE}}\n"
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
"{{OPENPMXTRAN_RECORD_FIELDS_DEFINE}}\n"
"\n"
"	/* declare and set IMODEL fields */\n"
"{{OPENPMXTRAN_IMODEL_FIELDS_DECLARE_AND_SET}}\n"
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
"{{OPENPMXTRAN_IMODEL_CODE}}\n"
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
"{{OPENPMXTRAN_IMODEL_FIELDS_SET}}\n"
"}\n"
"\n"
"typedef struct PREDICTVARS {\n"
"{{OPENPMXTRAN_PREDPARAMS_FIELDS_DECLARE}}\n"
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
"{{OPENPMXTRAN_RECORD_FIELDS_DEFINE}}\n"
"\n"
"	/* declare IMODEL fields */\n"
"{{OPENPMXTRAN_IMODEL_FIELDS_DEFINE}}\n"
"\n"
"#define THETA(i) 	((const double)_theta[(i)-1])\n"
"#define ETA(i) 	((const double)_eta[(i)-1])\n"
"#define STATE(i) 	((const double)_state[(i)-1])\n"
"#define A(i) 		((const double)_state[(i)-1])\n"
"#define DADT(i) 	(_dadt[(i)-1])\n"
"\n"
"{{OPENPMXTRAN_IMODEL_DIFFEQN_CODE}}\n"
"\n"
"#undef THETA\n"
"#undef ETA\n"
"#undef STATE\n"
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
"{{OPENPMXTRAN_RECORD_FIELDS_DEFINE}}\n"
"\n"
"	/* allow access to IMODEL fields */\n"
"{{OPENPMXTRAN_IMODEL_FIELDS_DEFINE}}\n"
"\n"
"	/* declare PREDICTVARS fields */\n"
"{{OPENPMXTRAN_PREDPARAMS_FIELDS_DECLARE}}\n"
"\n"
"#define THETA(i) 	((const double)_popparam->theta[(i)-1])\n"
"#define ETA(i) 	((const double)_popparam->eta[(i)-1])\n"
"#define ERR(i) 	((const double)_err[(i)-1])\n"
"#define EPS(i) 	((const double)_err[(i)-1])\n"
"#define A(i) 		((const double)_current->state[(i)-1])\n"
"\n"
"{{OPENPMXTRAN_PREDICT_CODE}}\n"
"\n"
"#undef THETA\n"
"#undef ETA\n"
"#undef A\n"
"#undef ERR\n"
"#undef EPS\n"
"\n"
"	/* set PREDICTVARS fields */\n"
"{{OPENPMXTRAN_PREDPARAMS_FIELDS_SET}}\n"
"\n"
"	return Y;\n"
"}\n"
"\n"
"static void openpmxtran_data_preprocess_callback(RECORD* _record)\n"
"{\n"
"	/* allow write access to RECORD fields */\n"
"{{OPENPMXTRAN_RECORD_FIELDS_WRITABLE_DEFINE}}\n"
"\n"
"	/* allow user to change data */\n"
"\tdo {\n"
"{{OPENPMXTRAN_DATA_PREPROCESS_CODE}}\n"
"\t} while (0);\n"
"\n"
"	/* write back RECORD fields */\n"
"{{OPENPMXTRAN_RECORD_FIELDS_WRITEBACK}}\n"
"\n"
"	return;\n"
"}\n"
"\n"
"static OPENPMX openpmx = (OPENPMX) {\n"
"	.filename = OPENPMXTRAN_DEFAULT_FILENAME,\n"
"	.nthread = 0,\n"
"	.data = { \n"
"{{OPENPMXTRAN_DATA_CONFIG}}\n"
"	},\n"
"	.advan = {\n"
"		.init = imodel_init,\n"
"		.predict = imodel_predict,\n"
"		.diffeqn = imodel_diffeqn,\n"
"		.imodelfields = {\n"
"			.size = sizeof(IMODEL),\n"
"			.field = {\n"
"{{OPENPMXTRAN_IMODEL_FIELDINFO}}\n"
"			},\n"
"		},\n"
"		.predictfields = {\n"
"			.size = sizeof(PREDICTVARS),\n"
"			.field = {\n"
"{{OPENPMXTRAN_PREDPARAMS_FIELDINFO}}\n"
"			},\n"
"		},\n"
"{{OPENPMXTRAN_ADVAN_INIT}}\n"
"	},\n"
"	.theta = {\n"
"{{OPENPMXTRAN_THETA_INIT}}\n"
"	},\n"
"	.omega = {\n"
"{{OPENPMXTRAN_OMEGA_INIT}}\n"
"	},\n"
"	.sigma = {\n"
"{{OPENPMXTRAN_SIGMA_INIT}}\n"
"	},\n"
"};\n"
"\n"
"OPENPMX current(void)\n"
"{\n"
"	return pmx_copy(&openpmx);\n"
"}\n"
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
"int main(void)\n"
"{\n"
"	/* do any preprocessing defined in the $DATA code */\n"
"	openpmxtran_data_preprocess();\n"
"//	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);\n"
"\n"
"{{OPENPMXTRAN_MAIN_CODE}}\n"
"\n"
"	pmx_cleanup(&openpmx);\n"
"	return EXIT_SUCCESS;\n"
"}\n";

