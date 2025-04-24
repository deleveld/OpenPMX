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
#include <string.h>

#include "openpmx.h"
#include "dataconfig.h"
#include "../utils/c22.h"

int structinfo_find_offset(const char* name, const STRUCTINFO* const structinfo)
{
	forcount(i, OPENPMX_FIELDS_MAX) {
		if (structinfo->field[i].name[0] == '\0')
			return -1;
		if (strcasecmp(name, structinfo->field[i].name) == 0)
			return structinfo->field[i].offset;
	}
	return -1;
}

