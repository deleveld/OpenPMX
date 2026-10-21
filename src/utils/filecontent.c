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

#include "filecontent.h"
#include "c22.h"

#include <stdio.h>
#include <stdlib.h>

void filecontent_free(FILECONTENT* file)
{
	free(file->mutptr);
}

FILECONTENT filecontent_alloc(const char* filename, ERRCTX* errctx)
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
	return (FILECONTENT) {
		.mutptr = buffer,
		.size = bytesRead,
	};

failed:
	if (buffer)
		free(buffer);
	if (fp)
		fclose(fp);
	return (FILECONTENT) { 0 };
}

