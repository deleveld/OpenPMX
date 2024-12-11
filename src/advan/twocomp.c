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
#include <math.h>

#include "advan.h"
#include "utils/c22.h"

typedef struct {
	ADVAN advan;
} ADVANCER_TWOCOMP;

static void advancer_twocomp_info(const struct ADVANFUNCS* const advanfuncs, FILE* f)
{
	(void) advanfuncs;

	fprintf(f, "advan model two compartment\n");
	fprintf(f, "advan parameters V1 V2 CL Q2\n");
}
static void advancer_twocomp_construct(ADVAN* advan, const ADVANFUNCS* const advanfuncs)
{
	ADVANCER_TWOCOMP* advantwocomp = (ADVANCER_TWOCOMP*)advan; /* cast up */

	assert(advanfuncs->advan_size == sizeof(ADVANCER_TWOCOMP));
	advan_base_construct(&advantwocomp->advan, advanfuncs); /* will zero entire object */
}

static void advancer_twocomp_destruct(ADVAN * advan)
{
	advan_base_destruct(advan);
}

typedef struct {
	ADVANFUNCS advanfuncs;

	int offsetV1;
	int offsetV2;
	int offsetCL;
	int offsetQ2;
} ADVANTABLE_TWOCOMP;

static void advancer_twocomp_advance_interval(ADVAN* advan,
											  const IMODEL* const imodel,
											  const RECORD* const record,
											  double* const state,
											  const POPPARAM* const popparam,
											  const double endtime,
											  const double* rates)
{
	/* no deed to cast up because we would immediately cast down again
	ADVANCER_TWOCOMP* advan = (ADVANCER_TWOCOMP*)advan; up cast */

	(void)record;
	(void)popparam;

	assert(advan->initcount > 0);

	const double time = advan->time;

	/* this advancer can only allow doses to first compartment */
	assert(endtime > time);
	assert(rates[0] >= 0.);
	assert(rates[1] == 0.);

	const ADVANTABLE_TWOCOMP* const imodeloffsets = (const ADVANTABLE_TWOCOMP*)advan->advanfuncs;
	const double V1 = *(const double*)(((char*)imodel) + imodeloffsets->offsetV1);
	const double V2 = *(const double*)(((char*)imodel) + imodeloffsets->offsetV2);
	const double CL = *(const double*)(((char*)imodel) + imodeloffsets->offsetCL);
	const double Q12 = *(const double*)(((char*)imodel) + imodeloffsets->offsetQ2);

	/* ADVAN 2-compartment code obtained from: Abuhelwa AY, Foster DJ, Upton RN. ADVAN-style analytical solutions for common pharmacokinetic imodels.
	   Journal of pharmacological and toxicological methods. 2015 Jun 30;73:42-8. */
	const double k10 = CL/V1;
	const double k12 = Q12/V1;
	const double k21 = k12*V1/V2;
	const double k20 = 0.;
	const double E1 = k10+k12;
	const double E2 = k21+k20;

    /* calculate hybrid rate constants */
    const double lambda1 = 0.5*((E1+E2)+sqrt(pow(E1+E2,2)-4.*(E1*E2-k12*k21)));
    const double lambda2 = 0.5*((E1+E2)-sqrt(pow(E1+E2,2)-4.*(E1*E2-k12*k21)));

#pragma push_macro("A1")
#pragma push_macro("A2")
#define A1 (state[0])
#define A2 (state[1])
	const double t = endtime - time; /* change in time */
	const double A1last = A1;
	const double A2last = A2;
	const double Doserate = rates[0];

	/* cached values for greater speed */
	const double exp_t_lambda1 = exp(-t*lambda1);
	const double exp_t_lambda2 = exp(-t*lambda2);

	const double A1term1 = (((A1last*E2+Doserate+A2last*k21)-A1last*lambda1)*exp_t_lambda1-((A1last*E2+Doserate+A2last*k21)-A1last*lambda2)*exp_t_lambda2)/(lambda2-lambda1);
	const double A1term2 = Doserate*E2*(1./(lambda1*lambda2)+exp_t_lambda1/(lambda1*(lambda1-lambda2))-exp_t_lambda2/(lambda2*(lambda1-lambda2)));
	A1 = A1term1+A1term2; /* Amount in the central compartment */

	const double A2term1 = (((A2last*E1+A1last*k12)-A2last*lambda1)*exp_t_lambda1-((A2last*E1+A1last*k12)-A2last*lambda2)*exp_t_lambda2)/(lambda2-lambda1);
	const double A2term2 = Doserate*k12*(1./(lambda1*lambda2)+exp_t_lambda1/(lambda1*(lambda1-lambda2))-exp_t_lambda2/(lambda2*(lambda1-lambda2)));
	A2 = A2term1+A2term2; /* Amount in the peripheral compartment */
#pragma pop_macro("A1")
#pragma pop_macro("A2")
}

ADVANFUNCS* pmx_advan_twocomp(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);
	assert(advanconfig->nstate == 0 || advanconfig->nstate == 2);

	let retinit = (ADVANTABLE_TWOCOMP) {
		.advanfuncs = {
			.advan_size = sizeof(ADVANCER_TWOCOMP),
			.construct = advancer_twocomp_construct,
			.destruct = advancer_twocomp_destruct,
			.info = advancer_twocomp_info,

			.reset = 0,
			.interval = advancer_twocomp_advance_interval,

			.advanconfig = advanconfig,
			.recordinfo = recordinfo_init(dataconfig),
			.nstate = 2,
		},
		.offsetV1 = structinfo_find_offset("V1", &advanconfig->imodelfields),
		.offsetV2 = structinfo_find_offset("V2", &advanconfig->imodelfields),
		.offsetCL = structinfo_find_offset("CL", &advanconfig->imodelfields),
		.offsetQ2 = structinfo_find_offset("Q2", &advanconfig->imodelfields),
	};
	assert(retinit.offsetV1 >= 0);
	assert(retinit.offsetV2 >= 0);
	assert(retinit.offsetCL >= 0);
	assert(retinit.offsetQ2 >= 0);

	ADVANTABLE_TWOCOMP* ret = malloc(sizeof(ADVANTABLE_TWOCOMP));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANTABLE_TWOCOMP));

	return &ret->advanfuncs;
}


