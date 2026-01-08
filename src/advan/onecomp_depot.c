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
 
/// This file implements an advancer for a one compartment model with a 
/// depot compartment.
///
/// The equations were obtained from: Abuhelwa AY, Foster DJ, Upton RN.
/// ADVAN-style analytical solutions for common pharmacokinetic imodels.
/// Journal of pharmacological and toxicological methods. 2015 Jun 30;73:42-8.

#include <assert.h>
#include <math.h>

#include "advan/advan.h"
#include "utils/c22.h"
#include "print.h"

typedef struct {
	ADVAN advan;
} ADVANCER_ONECOMP_DEPOT;

static void advancer_onecomp_depot_info(const struct ADVANFUNCS* const advanfuncs, FILE* f)
{
	(void) advanfuncs;

	fprintf(f, "advan model one compartment with absorbtion\n");
	fprintf(f, "advan parameters V CL KA\n");
}
static void advancer_onecomp_depot_construct(ADVAN* advan, const ADVANFUNCS* const advanfuncs)
{
	ADVANCER_ONECOMP_DEPOT* advanonecomp = (ADVANCER_ONECOMP_DEPOT*)advan; /* cast up */

	assert(advanfuncs->advan_size == sizeof(ADVANCER_ONECOMP_DEPOT));
	advan_base_construct(&advanonecomp->advan, advanfuncs); /* will zero entire object */
}

static void advancer_onecomp_depot_destruct(ADVAN * advan)
{
	advan_base_destruct(advan);
}

typedef struct {
	ADVANFUNCS advanfuncs;

	int offsetV;
	int offsetCL;
	int offsetKA;
} ADVANTABLE_ONECOMP_DEPOT;

static void advancer_onecomp_depot_advance_interval(ADVAN* advan,
													const IMODEL* const imodel,
													const RECORD* const record,
													double* const state,
													const POPPARAM* const popparam,
													const double endtime,
													const double* rates)
{
	/* no deed to cast up because we would immediately cast down again
	ADVANCER_ONECOMP_DEPOT* advan = (ADVANCER_ONECOMP_DEPOT*)advan; up cast */

	(void)record;
	(void)popparam;

	assert(advan->initcount > 0);

	const double time = advan->time;

	/* this advancer can only allow doses to first compartment */
	assert(endtime > time);
	assert(rates[0] >= 0.);
	assert(rates[1] == 0.);

	const ADVANTABLE_ONECOMP_DEPOT* const imodeloffsets = (const ADVANTABLE_ONECOMP_DEPOT*)advan->advanfuncs;
	const double V = *(const double*)(((char*)imodel) + imodeloffsets->offsetV);
	const double CL = *(const double*)(((char*)imodel) + imodeloffsets->offsetCL);
	const double KA = *(const double*)(((char*)imodel) + imodeloffsets->offsetKA);

	/* ADVAN 1-compartment code obtained from: Abuhelwa AY, Foster DJ, Upton RN. ADVAN-style analytical solutions for common pharmacokinetic imodels.
	   Journal of pharmacological and toxicological methods. 2015 Jun 30;73:42-8. */
	const double k10 = CL/V;

#pragma push_macro("A1")
#pragma push_macro("A2")
#define A1 (state[0])
#define A2 (state[1])
	const double t = endtime - time; /* change in time */
	double A1last = A1;
	double A2last = A2;
	const double Doserate = rates[0];
	
/// This advancer does not yet support infusions. See the GitHub issue. 
/// The workaround is to code the model as a differential equation.
	if (Doserate != 0.)
		fatal(0, "ADVANTABLE_ONECOMP_DEPOT fails with infusions, see GitHub issue. Use DIFFEQN.");

    A2last = A1last*KA/(KA-k10)*(exp(-t*k10)-exp(-t*KA))+A2last*exp(-t*k10);
    A1last = A1last*exp(-1.*t*KA);

    A2 = A2last;
    A1 = A1last;
    
#pragma pop_macro("A1")
#pragma pop_macro("A2")
}

ADVANFUNCS* pmx_advan_onecomp_depot(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);
	assert(advanconfig->nstate == 0 || advanconfig->nstate == 2);

	let retinit = (ADVANTABLE_ONECOMP_DEPOT) {
		.advanfuncs = {
			.advan_size = sizeof(ADVANCER_ONECOMP_DEPOT),
			.construct = advancer_onecomp_depot_construct,
			.destruct = advancer_onecomp_depot_destruct,
			.info = advancer_onecomp_depot_info,

			.reset = 0,
			.interval = advancer_onecomp_depot_advance_interval,

			.advanconfig = advanconfig,
			.recordinfo = recordinfo_init(dataconfig),
			.nstate = 2,
		},
		.offsetV = structinfo_find_offset("V", &advanconfig->imodelfields),
		.offsetCL = structinfo_find_offset("CL", &advanconfig->imodelfields),
		.offsetKA = structinfo_find_offset("KA", &advanconfig->imodelfields),
	};
	assert(retinit.offsetV >= 0);
	assert(retinit.offsetCL >= 0);
	assert(retinit.offsetKA >= 0);

	ADVANTABLE_ONECOMP_DEPOT* ret = malloc(sizeof(ADVANTABLE_ONECOMP_DEPOT));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANTABLE_ONECOMP_DEPOT));

	return &ret->advanfuncs;
}


