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

#include "getdelim.h"

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
    int c;
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



