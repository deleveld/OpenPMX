/*
 * This file is part of OpenPMX (https://github.com/deleveld/openpmx).
 * Copyright (c) 2022 Douglas Eleveld.
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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <values.h>
#include <unistd.h>
#include <stdarg.h>
#include <ctype.h>

#include "utils/c22.h"
#include "utils/vector.c"

typedef VECTOR(char*) STRINGS;

typedef VECTOR(char) STRING;

static void string_terminate(STRING* v)
{
	assert(v);
	vector_append(*v, 0);
	vector_remove_last(*v);
}

static void string_append(STRING* v, const char* t)
{
	let n = strlen(t);
	if (n > 0)
		vector_appendn(*v, t, n);
	string_terminate(v);
}

static void string_appendfv(STRING* v, const char* format, va_list args1)
{
	/* how much space do we need */
    va_list args2;
	va_copy(args2, args1);
    let nchars = vsnprintf(0, 0, format, args2);
    va_end(args2);

	/* reserve some space including zero terminator and do print at the old end */
	let oldsize = vector_size(*v);
    vector_resize(*v, oldsize + nchars + 1);
	char* oldend = &v->data[oldsize];
    vsnprintf(oldend, nchars + 1, format, args1);

    assert(vector_last(*v) == 0);	/* make sure we are terminated */
	vector_remove_last(*v);			/* remove terminator from vector */
}

static void string_appendf(STRING* v, const char* format, ... )
{
    va_list args1;

	va_start(args1, format);
	string_appendfv(v, format, args1);
	va_end(args1);
}

static STRING string_allocf(const char* format, ... )
{
    va_list args1;
	STRING ret = { 0 };

	va_start(args1, format);
	string_appendfv(&ret, format, args1);
	va_end(args1);

	return ret;
}

static void fatal(const char* format, ... )
{
    va_list args;

    va_start(args, format);
    fprintf(stderr, "error: ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    va_end(args);

    exit(EXIT_FAILURE);
}

static void strip_firstlast_space(char* s)
{
	let len = strlen(s);
	if (len == 0)
		return;

	/* remove space at end */
	char* a = s + len - 1;
	while (a != s && *a && isspace(*a)) {
		*a = 0;
		--a;
	}

	/* remove space at begin */
	a = s;
	while (*a && isspace(*a))
		++a;

	/* copy data over with classic strcpy
	 * which allows overlap */
	while ((*s++ = *a++)) { }
}

static int strip_quotes(char* s)
{
	let len = strlen(s);
	if (len == 0)
		return 1;

	/* remove quote at end */
	char* a = s + len - 1;
	if (*a != '\"')
		return 1;
	*a = 0;

	/* remove quote at begin */
	a = s;
	if (*a != '\"')
		return 1;
	++a;

	/* copy data over with classic strcpy
	 * which allows overlap */
	while ((*s++ = *a++)) { }

	return 0;
}

static void strip_comments(char *s, const char* ca, const char* cb, const char replace_char)
{
	let al = strlen(ca);
	let bl = strlen(cb);
	while (true) {
		char* a = strstr(s, ca);
		if (!a)
			break;
		char* b = strstr(a+al, cb);
		if (!b)
			break;
		b += bl;
		if (replace_char != 0)
			*a++ = replace_char;
		while ((*a++ = *b++)) { } /* classic strcpy */
	}
}

static char* read_file(const char* filename)
{
	var s = fopen(filename, "r");
	assert(s);

	/* https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c */
	char* buffer = NULL;
	size_t len;
	ssize_t bytes_read = getdelim(&buffer, &len, '\0', s);
	assert(bytes_read != -1);
	fclose(s);

	return buffer;
}

static int match_recordname(const char* p, const char* q, int* i, const bool whitedelim)
{
	var len = strlen(q);
	if (strncmp(p, q, len) == 0) {
		if (whitedelim && !isspace(p[len]))
			return 0;
		if (i)
			*i = len;
		return len;
	}
	return 0;
}

static void load_datafile_write_dataconfig(const char* datafile, STRING* data, STRINGS* fieldnames, STRING* config)
{
	fprintf(stdout, "data: \"%s\"\n", datafile);
	var stream = fopen(datafile, "r");
	if (!stream)
		fatal("could not open data");

	/* https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c */
	char* buffer = NULL;
	size_t len = 0;
	ssize_t bytes_read = getdelim(&buffer, &len, '\0', stream);
	fclose(stream);
	if (bytes_read == -1)
		fatal("could not read from data");

	/* allow comments in data files */
	strip_comments(buffer, "/*", "*/", ' ');
	strip_comments(buffer, "//", "\n", '\n');
	strip_firstlast_space(buffer);

	assert(vector_size(*fieldnames) == 0);

	string_append(data, "typedef struct RECORD {\n");
	string_append(data, "\tdouble ");

	/* process the file line by line */
	var nrows = -1;
	var line = buffer;
	var linebreak = strchr(line, '\n');
	if (!linebreak)
		fatal("could not any lines in data");
	*linebreak = 0;
	while (line) {
		const char delim[] = " \t,";

		/* read tokens in header */
		if (nrows == -1) {
			const char* token;
			char* rest = line;
			while ((token = strtok_r(rest, delim, &rest))) {

				/* write tokens to recordfields */
				let s = strdup(token);
				strip_firstlast_space(s);
				vector_append(*fieldnames, s);

				if (vector_size(*fieldnames) != 1)
					string_append(data, ", ");
				string_append(data, token);
			}
			string_append(data, ";\n");
			string_append(data, "} RECORD;\n");
			string_append(data, "RECORD OPENPMXTRAN_DATA_NAME[] = {\n");

		/* read tokens and write out in a RECORD struct */
		} else {
			var n = 0;
			char* token;
			char* rest = line;
			string_appendf(data, "\t{\t");
			while ((token = strtok_r(rest, delim, &rest))) {
				var s = strdup(token);
				strip_firstlast_space(s);
				if (strcmp(token, ".") == 0) {
					free(s);
					s = strdup("NAN");
				}
				string_append(data, s);
				free(s);
				if (n < vector_size(*fieldnames) - 1)
					string_append(data, ",");
				++n;
			}
			if (n != vector_size(*fieldnames))
				fatal("number fields in data line %i does not match header (%i)", n, vector_size(*fieldnames));
			string_append(data, "\t},\n");
		}
		++nrows;

		/* go to the next line, the last one wont have \n because we
		 * already stripped it */
		line = linebreak;
		if (line) {
			++line;
			linebreak = strchr(line, '\n');
			if (linebreak)
				*linebreak = 0;
		}
	}
	free(buffer);
	string_appendf(data, "};\n");

	/* Write out as an DATACONFIG initializer */
	string_append(config, "\t{\n");
	string_append(config, "\t\t.writeable = OPENPMXTRAN_DATA_NAME,\n");
	string_append(config, "\t\t.records = OPENPMXTRAN_DATA_NAME,\n");
	string_appendf(config, "\t\t.nrecords = %i,\n", nrows);
	string_appendf(config, "\t\t.offset_1 = true,\n");
	string_append(config, "\t\t.recordfields = {\n");
	string_append(config, "\t\t\t.size = sizeof(RECORD),\n");
	string_append(config, "\t\t\t.field = {\n");
	forvector(i, *fieldnames)
		string_appendf(config, "\t\t\t\t{ .name=\"%s\", .offset=offsetof(RECORD, %s) },\n", fieldnames->data[i], fieldnames->data[i]);
	string_append(config, "\t\t\t},\n");
	string_append(config, "\t\t},\n");
	string_append(config, "\t},");
}

typedef struct PARSERESULT {
	STRING filename;
	STRING data;
	STRING data_preprocess_code;
	STRING dataconfig;
	STRINGS record_field_names;

	STRING preamble;
	const char* advan_method;
	bool needs_diffeqn;
	const char* predict_state_maybe_unused;
	STRING advan_init;

	STRINGS imodel_field_names;
	STRING imodel_fields_code;

	STRING imodel_diffeqn_code;

	STRINGS predict_field_names;
	STRING predict_fields_code;

	STRING theta_init;
	STRING omega_init;
	STRING sigma_init;

	STRING analysis_code;

	STRING output;
} PARSERESULT;

static void parserresult_free(PARSERESULT* res)
{
	vector_free(res->filename);
	vector_free(res->data);
	vector_free(res->data_preprocess_code);
	vector_free(res->dataconfig);
	forvector(i, res->record_field_names)
		free(res->record_field_names.data[i]);
	vector_free(res->record_field_names);

	vector_free(res->preamble);
	vector_free(res->advan_init);

	vector_free(res->imodel_diffeqn_code);

	forvector(i, res->imodel_field_names)
		free(res->imodel_field_names.data[i]);
	vector_free(res->imodel_field_names);
	vector_free(res->imodel_fields_code);

	forvector(i, res->predict_field_names)
		free(res->predict_field_names.data[i]);
	vector_free(res->predict_field_names);
	vector_free(res->predict_fields_code);

	vector_free(res->theta_init);
	vector_free(res->omega_init);
	vector_free(res->sigma_init);

	vector_free(res->analysis_code);

	vector_free(res->output);
}

static char* match_brackets(char* p, int level)
{
	while (*p) {
		if (*p == '(')
			++level;
		else if (*p == ')') {
			--level;
			if (level == 0)
				return p;
		}
		++p;
	}
	return 0;
}

static void parse_data(PARSERESULT* res, char* p)
{
	if (vector_size(res->data) ||
		vector_size(res->record_field_names))
		fatal("more than one $DATA \"%s\"", p);

	var endvars = match_brackets(p, 1);
	if (!endvars)
		fatal("Could not find ')' to mark end of $DATA()");
	if (*endvars != ')')
		fatal("End of file within $DATA()");
	*endvars = 0;

	strip_firstlast_space(p);
	if (strip_quotes(p))
		fatal("$DATA filename missing quotes \"%s\"", p);

	load_datafile_write_dataconfig(p, &res->data, &res->record_field_names, &res->dataconfig);

	strip_firstlast_space(endvars + 1);
	string_append(&res->data_preprocess_code, endvars + 1);
}

static void parse_advan_init(PARSERESULT* res, char* p)
{
	if (vector_size(res->advan_init))
		fatal("more than one $ADVAN \"%s\"", p);

	var endvars = match_brackets(p, 1);
	if (!endvars)
		fatal("Could not find ')' to mark end of ADVAN() method name");
	*endvars = 0;

	strip_comments(p, "/*", "*/", ' ');
	strip_comments(p, "//", "\n", '\n');
	strip_firstlast_space(p);

	res->predict_state_maybe_unused = "";
	res->needs_diffeqn = false;
	if (streq(p, "pred")) {
		res->advan_method = "pmx_advan_pred";
		res->predict_state_maybe_unused = "	(void)_state;\n";
	} else if (streq(p, "onecomp")) 
		res->advan_method = "pmx_advan_onecomp";
	else if (streq(p, "onecomp_depot"))
		res->advan_method = "pmx_advan_onecomp_depot";
	else if (streq(p, "twocomp"))
		res->advan_method = "pmx_advan_twocomp";
	else if (streq(p, "threecomp"))
		res->advan_method = "pmx_advan_threecomp";
	else if (streq(p, "diffeqn_test")) {
		res->advan_method = "pmx_advan_diffeqn_test";
		res->needs_diffeqn = true;
	} else if (streq(p, "diffeqn_libgsl")) {
		res->advan_method = "pmx_advan_diffeqn_libgsl";
		res->needs_diffeqn = true;
	} else
		fatal("Could not identify ADVAN() method name");

	strip_firstlast_space(endvars + 1);
	string_append(&res->advan_init, endvars + 1);
}

static void parse_imodel(PARSERESULT* res, char* p)
{
	if (vector_size(res->imodel_field_names) ||
		vector_size(res->imodel_fields_code))
		fatal("more than one $IMODEL \"%s\"", p);

	/* parse field names */
	var endvars = match_brackets(p, 1);
	if (!endvars)
		fatal("could not find ')' to mark end of IMODEL() parameters");
	*endvars = 0;

	strip_comments(p, "/*", "*/", ' ');
	strip_comments(p, "//", "\n", '\n');
	strip_firstlast_space(p);

	char* token;
	char* rest = p;
	while ((token = strtok_r(rest, ", ", &rest))) {
		strip_firstlast_space(token);
		let s = strdup(token);
		vector_append(res->imodel_field_names, s);
	}

	/* after field names is code */
	strip_firstlast_space(endvars + 1);
	string_append(&res->imodel_fields_code, endvars + 1);
}

static void parse_predict(PARSERESULT* res, char* p)
{
	if (vector_size(res->predict_field_names) ||
		vector_size(res->predict_fields_code))
		fatal("more than one $PREDICT \"%s\"", p);

	/* parse field names */
	var endvars = match_brackets(p, 1);
	if (!endvars)
		fatal("could not find ')' to mark end of PREDICT() parameters");
	*endvars = 0;

	strip_comments(p, "/*", "*/", ' ');
	strip_comments(p, "//", "\n", '\n');
	strip_firstlast_space(p);

	char* token;
	char* rest = p;
	while ((token = strtok_r(rest, ", ", &rest))) {
		strip_firstlast_space(token);
		let s = strdup(token);
		vector_append(res->predict_field_names, s);
	}

	/* after field names is code */
	strip_firstlast_space(endvars + 1);
	string_append(&res->predict_fields_code, endvars + 1);
}

static void parse_diffeqn(PARSERESULT* res, char* p)
{
	if (vector_size(res->imodel_diffeqn_code))
		fatal("more than one $DIFFEQN \"%s\"", p);

	strip_firstlast_space(p);
	string_append(&res->imodel_diffeqn_code, p);
}

static void parse_theta(PARSERESULT* res, char* p)
{
	if (vector_size(res->theta_init))
		fatal("more than one $THETA \"%s\"", p);

	strip_firstlast_space(p);
	string_append(&res->theta_init, p);
}

static int count_elements(const char* p)
{
	var token = p;
	var n = 1;
	char* d;
	while ((d = strchr(token, ','))) {
		token = d + 1;
		++n;
	}
	return n;
}

static void parse_omega(PARSERESULT* res, char* p)
{
	var endvars = match_brackets(p, 1);
	if (!endvars)
		fatal("Could not find ')' to mark end of OMEGA()");
	*endvars = 0;

	strip_comments(p, "/*", "*/", ' ');
	strip_comments(p, "//", "\n", '\n');
	strip_firstlast_space(p);

	let n = count_elements(p);
	string_appendf(&res->omega_init, "\t\t{ OMEGA_DIAG, %i, { %s } },", n, p);
}

static void parse_omegablock(PARSERESULT* res, char* p)
{
	var endvars = match_brackets(p, 1);
	if (!endvars)
		fatal("Could not find ')' to mark end of OMEGABLOCK()");
	*endvars = 0;

	strip_comments(p, "/*", "*/", ' ');
	strip_comments(p, "//", "\n", '\n');
	strip_firstlast_space(p);

	/* calculate dimensions */
	/* https://www.wolframalpha.com/input?i2d=true&i=solve+k%5C%2844%29+n%3Dk*Divide%5B%5C%2840%29k%2B1%5C%2841%29%2C2%5D */
	let n = count_elements(p);
	let ndim = (int)floor((sqrt(8 * n + 1) - 1) / 2);

	/* TODO: maybe put some error checking on the dimension */
	string_appendf(&res->omega_init, "\t\t{ OMEGA_BLOCK, %i, { %s } },", ndim, p);
}

static void parse_omegasameblock(PARSERESULT* res, char* p)
{
	var endvars = match_brackets(p, 1);
	if (!endvars)
		fatal("Could not find ')' to mark end of OMEGASAME()");
	*endvars = 0;

	strip_comments(p, "/*", "*/", ' ');
	strip_comments(p, "//", "\n", '\n');
	strip_firstlast_space(p);

	string_appendf(&res->omega_init, "\t\t{ OMEGA_SAME, %s, { 0 } },", p);
}

static void parse_sigma(PARSERESULT* res, char* p)
{
	if (vector_size(res->sigma_init))
		fatal("more than one $SIGMA \"%s\"", p);

	var endvars = match_brackets(p, 1);
	if (!endvars)
		fatal("Could not find ')' to mark end of SIGMA()");
	if (*endvars != ')')
		fatal("End of file within $SIGMA()");
	*endvars = 0;

	strip_comments(p, "/*", "*/", ' ');
	strip_comments(p, "//", "\n", '\n');
	strip_firstlast_space(p);

	string_append(&res->sigma_init, p);

	var after = endvars + 1;
	strip_firstlast_space(after);
	if (*after != 0)
		fatal("Unknown text \"%s\" in $SIGMA", after);
}

static void parse_main(PARSERESULT* res, char* p)
{
	if (vector_size(res->analysis_code))
		fatal("More than one $MAIN \"%s\"", p);

	strip_firstlast_space(p);
	string_append(&res->analysis_code, p);
}

static int find_marker(const char* ptr)
{
	int i = 0;
	char prev_char = '\n';
	while (*ptr != 0) {
		if (*ptr == '$' && prev_char == '\n')
			break;

		/* look further */
		prev_char = *ptr;
		++ptr;
		++i;
	}
	return i;
}

int main(int argc, char* argv[])
{
	assert(argc == 2);
	const char* grfilename = argv[1];

	PARSERESULT res = { 0 };

	/* read in GR file */
	fprintf(stdout, "read: \"%s\"\n", grfilename);
	var grfile = read_file(grfilename);
	string_append(&res.filename, grfilename);

	/* everything up to first valid marker is a the preamble */
	var i = find_marker(grfile);
	char* begin = grfile + i;
	if (*begin == 0)
		fatal("could not find any markers in \"%s\"", grfilename);
	if (i > 0 && strlen(grfile)) {
		vector_appendn(res.preamble, grfile, i);
		string_terminate(&res.preamble);
	}

	while (*begin == '$') {
		var n = find_marker(begin + 1);
		char* end = begin + 1 + n;

		const char last_end = *end;
		*end = 0;

		int i;
		if (match_recordname(begin, "$DATA(", &i, false))
			parse_data(&res, begin + i);

		else if (match_recordname(begin, "$ADVAN(", &i, false))
			parse_advan_init(&res, begin + i);

		else if (match_recordname(begin, "$IMODEL(", &i, false))
			parse_imodel(&res, begin + i);

		else if (match_recordname(begin, "$PREDICT(", &i, false))
			parse_predict(&res, begin + i);

		else if (match_recordname(begin, "$DIFFEQN", &i, true))
			parse_diffeqn(&res, begin + i);

		else if (match_recordname(begin, "$THETA", &i, true))
			parse_theta(&res, begin + i);

		else if (match_recordname(begin, "$OMEGA(", &i, false))
			parse_omega(&res, begin + i);

		else if (match_recordname(begin, "$OMEGABLOCK(", &i, false))
			parse_omegablock(&res, begin + i);

		else if (match_recordname(begin, "$OMEGASAME(", &i, false))
			parse_omegasameblock(&res, begin + i);

		else if (match_recordname(begin, "$SIGMA(", &i, false))
			parse_sigma(&res, begin + i);

		else if (match_recordname(begin, "$MAIN", &i, true))
			parse_main(&res, begin + i);

		else
			fatal("unknown marker detected \"%s\"", begin);

		/* look for the next marker, stop if nothing to find */
		*end = last_end;
		begin = end;
	}
	free(grfile);

	/* copy over preamble */
	string_append(&res.output, "/* begin OPENPMXTRAN_PREAMBLE */\n");
	if (vector_size(res.preamble)) 
		string_append(&res.output, res.preamble.data);
	string_append(&res.output, "/* end OPENPMXTRAN_PREAMBLE */\n");

	/* read in template file */
extern char openpmxtran_template[];
	char* template_file = openpmxtran_template;

	/* everything up to first valid marker is a code */
	i = find_marker(template_file);
	begin = template_file + i;
	if (*begin == 0)
		fatal("could not find any markers in \"%s\"", grfilename);

	const char last_begin = *begin;
	*begin = 0;
	string_append(&res.output, template_file);
	*begin = last_begin;

	while (*begin == '$') {
		var n = find_marker(begin + 1);
		char* end = begin + 1 + n;

		const char last_end = *end;
		*end = 0;

		int i;
		if (match_recordname(begin, "${OPENPMXTRAN_DEFAULT_FILENAME}\n", &i, false)) {
			string_appendf(&res.output, "\"%s\"\n", res.filename.data);
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_DATA_ARRAY}\n", &i, false)) {
			string_appendf(&res.output, "%s\n", res.data.data);
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_DATA_CONFIG}\n", &i, false)) {
			string_appendf(&res.output, "%s\n", res.dataconfig.data);
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_IMODEL_FIELDS_DECLARE}\n", &i, false)) {
			forvector(i, res.imodel_field_names) {
				let s = res.imodel_field_names.data[i];
				string_appendf(&res.output, "\tdouble %s;\n", s);
			}
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_IMODEL_FIELDS_DECLARE_AND_SET}\n", &i, false)) {
			forvector(i, res.imodel_field_names) {
				let s = res.imodel_field_names.data[i];
				string_appendf(&res.output, "\tdouble %s = _imodel->%s;\n", s, s);
			}
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_RECORD_FIELDS_DEFINE}\n", &i, false)) {
			forvector(i, res.record_field_names) {
				let s = res.record_field_names.data[i];
				string_appendf(&res.output, "\tconst double %s = _record->%s;\n", s, s);
				string_appendf(&res.output, "\t(void)%s;\n", s);
			}
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_RECORD_FIELDS_WRITABLE_DEFINE}\n", &i, false)) {
			forvector(i, res.record_field_names) {
				let s = res.record_field_names.data[i];
				string_appendf(&res.output, "\tdouble %s = _record->%s;\n", s, s);
			}
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_DATA_PREPROCESS_CODE}\n", &i, false)) {
			if (vector_size(res.data_preprocess_code))
				string_appendf(&res.output, "%s\n", res.data_preprocess_code.data);
			else
				string_appendf(&res.output, "	/* no data preprocess */\n");

			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_RECORD_FIELDS_WRITEBACK}\n", &i, false)) {
			forvector(i, res.record_field_names) {
				let s = res.record_field_names.data[i];
				string_appendf(&res.output, "\t_record->%s = %s;\n", s, s);
			}
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_IMODEL_FIELDS_DECLARE}\n", &i, false)) {
			forvector(i, res.imodel_field_names) {
				let s = res.imodel_field_names.data[i];
				string_appendf(&res.output, "\tdouble %s;\n", s);
			}
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_IMODEL_FIELDS_CODE}\n", &i, false)) {
			string_appendf(&res.output, "\t%s\n", res.imodel_fields_code.data);
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_IMODEL_FIELDS_SET}\n", &i, false)) {
			forvector(i, res.imodel_field_names) {
				let s = res.imodel_field_names.data[i];
				string_appendf(&res.output, "\t_imodel-> %s = %s;\n", s, s);
			}
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_IMODEL_FIELDS_DECLARE}\n", &i, false)) {
			forvector(i, res.imodel_field_names) {
				let s = res.imodel_field_names.data[i];
				string_appendf(&res.output, "\tdouble %s;\n", s);
			}
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_IMODEL_FIELDS_DEFINE}\n", &i, false)) {
			forvector(i, res.imodel_field_names) {
				let s = res.imodel_field_names.data[i];
				string_appendf(&res.output, "\tconst double %s = _imodel->%s;\n", s, s);
				string_appendf(&res.output, "\t(void)%s;\n", s);
			}
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_PREDPARAMS_FIELDS_DECLARE}\n", &i, false)) {
			forvector(i, res.predict_field_names) {
				let s = res.predict_field_names.data[i];
				string_appendf(&res.output, "\tdouble %s;\n", s);
			}
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_IMODEL_DIFFEQN_CODE}\n", &i, false)) {
			if (vector_size(res.imodel_diffeqn_code)) {
				if (res.needs_diffeqn == false)
					fatal("This $ADVAN does not require $DIFFEQN");
				string_appendf(&res.output, "\t%s\n", res.imodel_diffeqn_code.data);
			} else
				string_appendf(&res.output, "\t(void)_dadt;\n\t/* no $DIFFEQN code*/\n\texit(EXIT_FAILURE);\n");

			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_PREDICT_CODE}\n", &i, false)) {
			string_append(&res.output, res.predict_state_maybe_unused);
			string_appendf(&res.output, "\t%s\n", res.predict_fields_code.data);
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_PREDPARAMS_FIELDS_SET}\n", &i, false)) {
			forvector(i, res.predict_field_names) {
				let s = res.predict_field_names.data[i];
				string_appendf(&res.output, "\t_predparams-> %s = %s;\n", s, s);
			}
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_ADVAN_INIT}\n", &i, false)) {
			string_appendf(&res.output, "\t\t.method = %s,\n", res.advan_method);
			string_appendf(&res.output, "\t\t%s\n", res.advan_init.data);
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_IMODEL_FIELDINFO}\n", &i, false)) {
			forvector(i, res.imodel_field_names) {
				let s = res.imodel_field_names.data[i];
				string_appendf(&res.output, "\t\t\t\t{ .name=\"%s\", .offset = offsetof(IMODEL, %s) },\n", s, s);
			}
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_PREDPARAMS_FIELDINFO}\n", &i, false)) {
			forvector(i, res.predict_field_names) {
				let s = res.predict_field_names.data[i];
				string_appendf(&res.output, "\t\t\t\t{ .name=\"%s\", .offset = offsetof(PREDICTVARS, %s) },\n", s, s);
			}
			string_append(&res.output, begin + i);

		} else if (match_recordname(begin, "${OPENPMXTRAN_THETA_INIT}\n", &i, false)) {
			if (vector_size(res.theta_init)) 
				string_appendf(&res.output, "\t%s\n", res.theta_init.data);
			string_append(&res.output, begin + i);
		} else if (match_recordname(begin, "${OPENPMXTRAN_OMEGA_INIT}\n", &i, false)) {
			if (vector_size(res.omega_init)) 
				string_appendf(&res.output, "\t%s\n", res.omega_init.data);
			string_append(&res.output, begin + i);
		} else if (match_recordname(begin, "${OPENPMXTRAN_SIGMA_INIT}\n", &i, false)) {
			if (vector_size(res.sigma_init)) 
				string_appendf(&res.output, "\t%s\n", res.sigma_init.data);
			string_append(&res.output, begin + i);
		
		} else if (match_recordname(begin, "${OPENPMXTRAN_MAIN_CODE}\n", &i, false)) {
			string_appendf(&res.output, "\t%s\n", res.analysis_code.data);
			string_append(&res.output, begin + i);

		} else
			string_append(&res.analysis_code, begin);

		/* look for the next marker, stop if nothing to find */
		*end = last_end;
		begin = end;
	}

	/* write output file */
	STRING outputfile = string_allocf("%s.c", res.filename.data);
	fprintf(stdout, "write: \"%s\"\n", outputfile.data);
	var f = fopen(outputfile.data, "w");
	assert(f);
	fputs(res.output.data, f);
	fclose(f);
	vector_free(outputfile);

	parserresult_free(&res);

	return EXIT_SUCCESS;
}

char openpmxtran_template[] =
"/*\n"
" * Code below this message is generated by openpmx (https://github.com/deleveld/openpmx).\n"
" * Copyright (c) 2022 Douglas Eleveld.\n"
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
"#include \"openpmx.h\"\n"
"\n"
"#define OPENPMXTRAN_DEFAULT_FILENAME \\\n"
"${OPENPMXTRAN_DEFAULT_FILENAME}\n"
"\n"
"/* begin OPENPMXTRAN_DATA_ARRAY */\n"
"${OPENPMXTRAN_DATA_ARRAY}\n"
"/* end OPENPMXTRAN_DATA_ARRAY */\n"
"\n"
"typedef struct IMODEL {\n"
"/* begin OPENPMXTRAN_IMODEL_FIELDS_DECLARE */\n"
"${OPENPMXTRAN_IMODEL_FIELDS_DECLARE}\n"
"/* end OPENPMXTRAN_IMODEL_FIELDS_DECLARE */\n"
"} IMODEL;\n"
"\n"
"static void imodel_init(IMODEL* const _imodel,\n"
"						 ADVANSTATE* const _advanstate,\n"
"						 const POPPARAM* const _popparam)\n"
"{\n"
"	const double* _theta = _popparam->theta;\n"
"	const double* _eta = _popparam->eta;\n"
"	const RECORD* const _record = _advanstate->record;\n"
"	const int initcount = _advanstate->initcount;\n"
"	(void) initcount;\n"
"\n"
"	/* allow access to RECORD fields */\n"
"${OPENPMXTRAN_RECORD_FIELDS_DEFINE}\n"
"\n"
"	/* declare and set IMODEL fields */\n"
"${OPENPMXTRAN_IMODEL_FIELDS_DECLARE_AND_SET}\n"
"\n"
"#define THETA(i) 		((const double)_theta[i-1])\n"
"#define ETA(i) 		((const double)_eta[i-1])\n"
"#define A(i) 			(_advanstate->state[i-1])\n"
"#define ALAG(i) 		(_advanstate->amtlag[i-1])\n"
"#define F(i) 			(_advanstate->bioavail[i-1])\n"
"#define STATETIME 		(_advanstate->statetime)\n"
"#define INITTIME(t) 	pmx_advan_inittime(_advanstate, (t))\n"
"\n"
"/* begin user code to init the IMODEL */\n"
"${OPENPMXTRAN_IMODEL_FIELDS_CODE}\n"
"/* end user code to init the IMODEL */\n"
"\n"
"#undef THETA\n"
"#undef ETA\n"
"#undef STATE\n"
"#undef A\n"
"#undef ALAG\n"
"#undef STATETIME\n"
"#undef INITTIME\n"
"#undef F\n"
"\n"
"	/* set IMODEL fields */\n"
"${OPENPMXTRAN_IMODEL_FIELDS_SET}\n"
"}\n"
"\n"
"typedef struct PREDICTVARS {\n"
"${OPENPMXTRAN_PREDPARAMS_FIELDS_DECLARE}\n"
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
"${OPENPMXTRAN_RECORD_FIELDS_DEFINE}\n"
"\n"
"	/* declare IMODEL fields */\n"
"${OPENPMXTRAN_IMODEL_FIELDS_DEFINE}\n"
"\n"
"#define THETA(i) 	((const double)_theta[i-1])\n"
"#define ETA(i) 	((const double)_eta[i-1])\n"
"#define STATE(i) 	((const double)_state[i-1])\n"
"#define A(i) 		((const double)_state[i-1])\n"
"#define DADT(i) 	(_dadt[i-1])\n"
"\n"
"/* begin user code to do the differential equations */\n"
"${OPENPMXTRAN_IMODEL_DIFFEQN_CODE}\n"
"/* end user code to do the differential equations */\n"
"\n"
"#undef THETA\n"
"#undef ETA\n"
"#undef STATE\n"
"#undef A\n"
"#undef DADT\n"
"}\n"
"\n"
"static double imodel_predict(const IMODEL* const _imodel,\n"
"							 const PREDICTSTATE* const _predictstate,\n"
"							 const POPPARAM* const _popparam,\n"
"							 const double* const _err,\n"
"							 PREDICTVARS* _predparams)\n"
"{\n"
"	const RECORD* const _record = _predictstate->record;\n"
"	const double* const _state = _predictstate->state;\n"
"	(void) _popparam;\n"
"	(void) _record;\n"
"	(void) _predparams;\n"
"	double Y = NAN;\n"
"\n"
"	/* allow access to RECORD fields */\n"
"${OPENPMXTRAN_RECORD_FIELDS_DEFINE}\n"
"\n"
"	/* allow access to IMODEL fields */\n"
"${OPENPMXTRAN_IMODEL_FIELDS_DEFINE}\n"
"\n"
"	/* declare PREDICTVARS fields */\n"
"${OPENPMXTRAN_PREDPARAMS_FIELDS_DECLARE}\n"
"\n"
"#define THETA(i) 	((const double)_theta[i-1])\n"
"#define ETA(i) 	((const double)_popparam->eta[i-1])\n"
"#define A(i) 		((const double)_state[i-1])\n"
"#define ERR(i) 	((const double)_err[i-1])\n"
"#define EPS(i) 	((const double)_err[i-1])\n"
"\n"
"/* begin user code to do prediction, set Y and PREDICTVARS */\n"
"${OPENPMXTRAN_PREDICT_CODE}\n"
"/* end user code to do prediction, set Y and PREDICTVARS */\n"
"\n"
"#undef THETA\n"
"#undef ETA\n"
"#undef STATE\n"
"#undef A\n"
"#undef ERR\n"
"\n"
"	/* set PREDICTVARS fields */\n"
"${OPENPMXTRAN_PREDPARAMS_FIELDS_SET}\n"
"\n"
"	return Y;\n"
"}\n"
"\n"
"static bool openpmxtran_data_preprocess_callback(RECORD* _record)\n"
"{\n"
"	/* set to remove this record from dataset analysis */\n"
"	bool remove = false;\n"
"\n"
"	/* allow write access to RECORD fields */\n"
"${OPENPMXTRAN_RECORD_FIELDS_WRITABLE_DEFINE}\n"
"\n"
"	/* allow user to change data */\n"
"${OPENPMXTRAN_DATA_PREPROCESS_CODE}\n"
"\n"
"	/* write back RECORD fields */\n"
"${OPENPMXTRAN_RECORD_FIELDS_WRITEBACK}\n"
"\n"
"	return remove ? false : true;\n"
"}\n"
"\n"
"static OPENPMX openpmx = (OPENPMX) {\n"
"	.filename = OPENPMXTRAN_DEFAULT_FILENAME,\n"
"	.nthread = 0,\n"
"	.offset_1 = true,\n"
"	.data = \n"
"/* begin OPENPMXTRAN_DATA_CONFIG */\n"
"${OPENPMXTRAN_DATA_CONFIG}\n"
"/* end OPENPMXTRAN_DATA_CONFIG */\n"
"	.advan = {\n"
"		.init = imodel_init,\n"
"		.predict = imodel_predict,\n"
"		.diffeqn = imodel_diffeqn,\n"
"		.imodelfields = {\n"
"			.size = sizeof(IMODEL),\n"
"			.field = {\n"
"/* begin OPENPMXTRAN_IMODEL_FIELDINFO */\n"
"${OPENPMXTRAN_IMODEL_FIELDINFO}\n"
"/* end OPENPMXTRAN_IMODEL_FIELDINFO */\n"
"			},\n"
"		},\n"
"		.predictfields = {\n"
"			.size = sizeof(PREDICTVARS),\n"
"			.field = {\n"
"/* begin OPENPMXTRAN_PREDPARAMS_FIELDINFO */\n"
"${OPENPMXTRAN_PREDPARAMS_FIELDINFO}\n"
"/* end OPENPMXTRAN_PREDPARAMS_FIELDINFO */\n"
"			},\n"
"		},\n"
"/* begin OPENPMXTRAN_ADVAN_INIT */\n"
"${OPENPMXTRAN_ADVAN_INIT}\n"
"/* end OPENPMXTRAN_ADVAN_INIT */\n"
"	},\n"
"	.theta = {\n"
"/* begin OPENPMXTRAN_THETA_INIT */\n"
"${OPENPMXTRAN_THETA_INIT}\n"
"/* end OPENPMXTRAN_THETA_INIT */\n"
"	},\n"
"	.omega = {\n"
"/* begin OPENPMXTRAN_OMEGA_INIT */\n"
"${OPENPMXTRAN_OMEGA_INIT}\n"
"/* end OPENPMXTRAN_OMEGA_INIT */\n"
"	},\n"
"	.sigma = {\n"
"/* begin OPENPMXTRAN_SIGMA_INIT */\n"
"${OPENPMXTRAN_SIGMA_INIT}\n"
"/* end OPENPMXTRAN_SIGMA_INIT */\n"
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
"static void openpmxtran_data_preprocess(void)\n"
"{\n"
"	RECORD* data = openpmx.data.writeable;\n"
"	const int nrecords = openpmx.data.nrecords;\n"
"	int i = 0;\n"
"	int n = 0;\n"
"	for (i=0; i<nrecords; i++) {\n"
"		RECORD r = data[i];\n"
"		if (openpmxtran_data_preprocess_callback(&r)) {\n"
"			data[n] = r;\n"
"			++n;\n"
"		}\n"
"	}\n"
"	for (i=n; i<nrecords; i++)\n"
"		data[n].ID = NAN;\n"
"}\n"
"\n"
"int main(void)\n"
"{\n"
"	/* do any preprocessing defined in the $DATA code */\n"
"	openpmxtran_data_preprocess();\n"
"\n"
"/* begin OPENPMXTRAN_MAIN_CODE */\n"
"${OPENPMXTRAN_MAIN_CODE}\n"
"/* end OPENPMXTRAN_MAIN_CODE */\n"
"\n"
"	pmx_cleanup(&openpmx);\n"
"	return EXIT_SUCCESS;\n"
"}\n";
