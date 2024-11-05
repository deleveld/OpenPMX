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

#include <assert.h>
#include <math.h>

#include "advan.h"
#include "utils/c22.h"

typedef struct {
	ADVAN advan;
} ADVANCER_ONECOMP;

static void advancer_onecomp_info(const struct ADVANFUNCS* const advanfuncs, FILE* f)
{
	(void) advanfuncs;

	fprintf(f, "advan model: one compartment\n");
	fprintf(f, "advan parameters: V1 CL\n");
}
static void advancer_onecomp_construct(ADVAN* advan, const ADVANFUNCS* const advanfuncs)
{
	ADVANCER_ONECOMP* advanonecomp = (ADVANCER_ONECOMP*)advan; /* cast up */

	assert(advanfuncs->advan_size == sizeof(ADVANCER_ONECOMP));
	advan_base_construct(&advanonecomp->advan, advanfuncs); /* will zero entire object */
}

static void advancer_onecomp_destruct(ADVAN * advan)
{
	advan_base_destruct(advan);
}

typedef struct {
	ADVANFUNCS advanfuncs;

	int offsetV1;
	int offsetCL;
} ADVANTABLE_ONECOMP;

static void advancer_onecomp_advance_interval(ADVAN* advan,
											  const IMODEL* const imodel,
											  const RECORD* const record,
											  double* const state,
											  const POPPARAM* const popparam,
											  const double endtime,
											  const double* rates)
{
	/* no deed to cast up because we would immediately cast down again
	ADVANCER_ONECOMP* advan = (ADVANCER_ONECOMP*)advan; up cast */

	(void)record;
	(void)popparam;

	assert(advan->initcount > 0);

	const double time = advan->time;

	/* this advancer can only allow doses to first compartment */
	assert(endtime > time);
	assert(rates[0] >= 0.);
	assert(rates[1] == 0.);

	const ADVANTABLE_ONECOMP* const imodeloffsets = (const ADVANTABLE_ONECOMP*)advan->advanfuncs;
	const double V1 = *(const double*)(((char*)imodel) + imodeloffsets->offsetV1);
	const double CL = *(const double*)(((char*)imodel) + imodeloffsets->offsetCL);

	/* ADVAN 1-compartment code obtained from: Abuhelwa AY, Foster DJ, Upton RN. ADVAN-style analytical solutions for common pharmacokinetic imodels.
	   Journal of pharmacological and toxicological methods. 2015 Jun 30;73:42-8. */
	const double k10 = CL/V1;

#pragma push_macro("A1")
#define A1 (state[0])
	const double t = endtime - time; /* change in time */
	const double A1last = A1;
	const double Doserate = rates[0];

    A1 = Doserate/k10*(1-exp(-t*k10))+A1last*exp(-t*k10);
#pragma pop_macro("A1")
}

ADVANFUNCS* pmx_advan_onecomp(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);
	assert(advanconfig->nstate == 0 || advanconfig->nstate == 1);

	let retinit = (ADVANTABLE_ONECOMP) {
		.advanfuncs = {
			.advan_size = sizeof(ADVANCER_ONECOMP),
			.construct = advancer_onecomp_construct,
			.destruct = advancer_onecomp_destruct,
			.info = advancer_onecomp_info,

			.reset = 0,
			.interval = advancer_onecomp_advance_interval,

			.advanconfig = advanconfig,
			.recordinfo = recordinfo_init(dataconfig),
			.nstate = 1,
		},
		.offsetV1 = structinfo_find_offset("V1", &advanconfig->imodelfields),
		.offsetCL = structinfo_find_offset("CL", &advanconfig->imodelfields),
	};
	assert(retinit.offsetV1 >= 0);
	assert(retinit.offsetCL >= 0);

	ADVANTABLE_ONECOMP* ret = malloc(sizeof(ADVANTABLE_ONECOMP));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANTABLE_ONECOMP));

	return &ret->advanfuncs;
}


