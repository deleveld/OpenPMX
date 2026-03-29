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

/// This file define the file extenstions and numeric formats used.

#ifndef OPENPMX_DEFINES_H
#define OPENPMX_DEFINES_H

#ifdef __cplusplus
extern "C" {
#endif

#define OPENPMX_OUTFILE 				".out"
#define OPENPMX_EXTFILE 				".ext"
#define OPENPMX_PHIFILE 				".phi"
#define OPENPMX_YHATFILE 				".yhat"
#define OPENPMX_ICOVRESAMPLEFILE 		".icov"
#define OPENPMX_TABLEFILE				".txt"

#define OPENPMX_SFORMAT 				" %12s"
#define OPENPMX_IFORMAT 				" %12i"
#define OPENPMX_FFORMAT 				" %12g"

#define OPENPMX_HEADER_FORMAT 			" %24s"
#define OPENPMX_TABLE_FORMAT 			" %24.17g"

#define OPENPMX_IMODEL_MAX				64
#define OPENPMX_PREDICTVARS_MAX			64
#define OPENPMX_STATE_MAX				32
#define OPENPMX_SIMULINFUSION_MAX		16

#define OPENPMX_EXTFILE_LOWER_BOUNDS	-2000000001
#define OPENPMX_EXTFILE_UPPER_BOUNDS	-2000000002
#define OPENPMX_EXTFILE_EXTRA_FIXED		-2000000006
#define OPENPMX_EXTFILE_OMEGA_BLOCKS	-2000000007

#define OPENPMX_SERVER_FLUSH_PERIOD_MS	5000
#define OPENPMX_SERVER_TIMEOUT_S		(15*60)

#ifdef __cplusplus
}
#endif

#endif
