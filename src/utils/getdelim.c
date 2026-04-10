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
#include <ctype.h>

#include "getdelim.h"
#include "c22.h"

ssize_t openpmx_getdelim(char **lineptr, size_t *n, int delim, FILE *stream)
{
    if (lineptr == NULL || n == NULL || stream == NULL)
        return -1;

    if (*lineptr == NULL || *n == 0) {
        *n = 512;  // initial size
        *lineptr = malloc(*n);
        if (*lineptr == NULL)
            return -1;
    }

    size_t pos = 0;
    int c = EOF;
    while ((c = fgetc(stream)) != EOF) {
        if (pos + 1 >= *n) {
            size_t new_size = *n * 2;
            char *new_ptr = realloc(*lineptr, new_size);
            if (!new_ptr)
                return -1;
            *lineptr = new_ptr;
            *n = new_size;
        }
        (*lineptr)[pos++] = (char)c;
        if (c == delim)
            break;
    }
    if (pos == 0 && c == EOF)  // nothing read
        return -1;
    (*lineptr)[pos] = '\0';
    return pos;
}

void strip_firstlast_space(char* s)
{
	/* remove space at begin */
	char* begin = s;
	while (*begin && isspace((unsigned char)*begin))
		++begin;

	/* remove space at end */
	size_t len = strlen(begin);
	while (len > 0 && isspace((unsigned char)begin[len - 1]))
		--len;

	/* copy over, must allow overlap */
	memmove(s, begin, len);
	s[len] = '\0';
}

static int get_delim_tokens_comma(char* line, VECPTR* namevec)
{
    if (line[0] == ',')
        return 1;

    char *start = line;
    char *end;

    while ((end = strchr(start, ',')) != NULL) {
        *end = '\0';
        strip_firstlast_space(start);
        if (strlen(start)) {
            vector_append(*namevec, start);
        } else
            return 1;   /* missing token */
        start = end + 1;
    }

    /* Handle the last (or only) token after the final comma */
    strip_firstlast_space(start);
    if (strlen(start)) {
        vector_append(*namevec, start);
    } else
        return 1;   /* missing token */

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

int get_delim_tokens(char* line, VECPTR* namevec, const GET_DELIM_SEP sep)
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
		return 2;
	}
		
	return 0;
}


