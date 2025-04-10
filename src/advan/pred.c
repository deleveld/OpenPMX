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

#include <assert.h>
#include <string.h>

#include "advan.h"
#include "utils/c22.h"

typedef struct {
	ADVAN advancer_object;
} ADVANCER_PRED;

static void advancer_pred_info(const struct ADVANFUNCS* const advanfuncs, FILE* f)
{
	(void) advanfuncs;

	fprintf(f, "advan model pred\n");
}
static void advancer_pred_construct(ADVAN* advan, const struct ADVANFUNCS* const advanfuncs)
{
	/* no need to cast up since the only thing we do is cast down again
	ADVANCER_PRED* advan = (ADVANCER_PRED*)advan; cast up */

	assert(advanfuncs->advan_size == sizeof(ADVANCER_PRED));

	advan_base_construct(advan, advanfuncs);
}

static void advancer_pred_destruct(ADVAN * advan)
{
	advan_base_destruct(advan);
}

static void advancer_pred_advance_interval(ADVAN* advan,
										   const IMODEL* const imodel,
										   const RECORD* const record,
										   double* const state,
										   const POPPARAM* const popparam,
										   const double endtime,
										   const double* rates)
{
	(void)imodel;
	(void)record;
	(void)state;
	(void)popparam;
	(void)rates;

	advan->time = endtime;
}

ADVANFUNCS* pmx_advan_pred(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);
	assert(advanconfig->nstate == 0);

	let retinit = (ADVANFUNCS) {
		.advan_size = sizeof(ADVANCER_PRED),
		.construct = advancer_pred_construct,
		.destruct = advancer_pred_destruct,
		.info = advancer_pred_info,

		.reset = 0,
		.interval = advancer_pred_advance_interval,

		.advanconfig = advanconfig,
		.recordinfo = recordinfo_init(dataconfig),
		.nstate = 0,
	};

	/* make binary copy so init can have const members */
	ADVANFUNCS* ret = malloc(sizeof(ADVANFUNCS));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANFUNCS));

	return ret;
}
