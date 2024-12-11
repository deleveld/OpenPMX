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

/* ADVAN 3-compartment code obtained from: Abuhelwa AY, Foster DJ, Upton RN.
 * ADVAN-style analytical solutions for common pharmacokinetic imodels.
   Journal of pharmacological and toxicological methods. 2015 Jun 30;73:42-8. */

#include <assert.h>
#include <math.h>

#include "advan.h"
#include "utils/c22.h"

typedef struct {
	ADVAN advan;

	int last_reset_count;

	double k12;
	double k21;
	double k13;
	double k31;
	double E1;
	double E2;
	double E3;
	double lambda1;
	double lambda2;
	double lambda3;
} ADVANCER_THREECOMP;

static void advancer_threecomp_info(const struct ADVANFUNCS* const advanfuncs,
									FILE* f)
{
	(void) advanfuncs;

	fprintf(f, "advan model three compartment\n");
	fprintf(f, "advan parameters V1 V2 V3 CL Q2 Q3\n");
}

void advancer_threecomp_construct(ADVAN* advan,
								  const struct ADVANFUNCS* const advanfuncs)
{
	ADVANCER_THREECOMP* advanthreecomp = (ADVANCER_THREECOMP*)advan; /* cast up */

	assert(advanfuncs->advan_size == sizeof(ADVANCER_THREECOMP));
	advan_base_construct(&advanthreecomp->advan, advanfuncs); /* will zero full size */

	advanthreecomp->last_reset_count = 0;
}

static void advancer_threecomp_destruct(ADVAN * advan)
{
	advan_base_destruct(advan);
}

typedef struct {
	ADVANFUNCS advanfuncs;

	int offsetV1;
	int offsetV2;
	int offsetV3;
	int offsetCL;
	int offsetQ2;
	int offsetQ3;
} ADVANTABLE_THREECOMP;

static inline void imodel_threecomp_reset_macro_constants(ADVANCER_THREECOMP* advanthreecomp, const IMODEL* imodel)
{
	const ADVANTABLE_THREECOMP* const imodeloffsets = (const ADVANTABLE_THREECOMP*)advanthreecomp->advan.advanfuncs;

	const double V1 = *(const double*)(((char*)imodel) + imodeloffsets->offsetV1);
	const double V2 = *(const double*)(((char*)imodel) + imodeloffsets->offsetV2);
	const double V3 = *(const double*)(((char*)imodel) + imodeloffsets->offsetV3);
	const double CL = *(const double*)(((char*)imodel) + imodeloffsets->offsetCL);
	const double Q12 = *(const double*)(((char*)imodel) + imodeloffsets->offsetQ2);
	const double Q13 = *(const double*)(((char*)imodel) + imodeloffsets->offsetQ3);
	const double k10 = CL/V1;
	const double k12 = Q12/V1;
	const double k21 = k12*V1/V2;
	const double k13 = Q13/V1;
	const double k31 = k13*V1/V3;
	const double k20 = 0.;
	const double k30 = 0.;
	const double E1 = k10+k12+k13;
	const double E2 = k21+k20;
	const double E3 = k31+k30;

	/* calculate hybrid rate constants */
	const double a = E1+E2+E3;
	const double b = E1*E2+E3*(E1+E2)-k12*k21-k13*k31;
	const double c = E1*E2*E3-E3*k12*k21-E2*k13*k31;

	const double m = (3.*b - pow(a,2.))/3.;
	const double n = (2.*pow(a,3.) - 9.*a*b + 27.*c)/27.;
	const double Q = (pow(n,2.))/4. + (pow(m,3.))/27.;

	const double alpha = sqrt(-1.*Q);
	const double beta = -1.*n/2.;
	const double gamma = sqrt(pow(beta,2.)+pow(alpha,2.));
	const double theta = atan2(alpha,beta);

	const double pow_gamma_1_3 = cbrt(gamma);
	const double lambda1 = a/3. + pow_gamma_1_3*(cos(theta/3.) + sqrt(3.)*sin(theta/3.));
	const double lambda2 = a/3. + pow_gamma_1_3*(cos(theta/3.) - sqrt(3.)*sin(theta/3.));
	const double lambda3 = a/3. -(2.*pow_gamma_1_3*cos(theta/3.));

	advanthreecomp->k12 = k12;
	advanthreecomp->k21 = k21;
	advanthreecomp->k13 = k13;
	advanthreecomp->k31 = k31;
	advanthreecomp->E1 = E1;
	advanthreecomp->E2 = E2;
	advanthreecomp->E3 = E3;
	advanthreecomp->lambda1 = lambda1;
	advanthreecomp->lambda2 = lambda2;
	advanthreecomp->lambda3 = lambda3;
}

static void advancer_threecomp_advance_interval(ADVAN* advan,
												const IMODEL* const imodel,
												const RECORD* const record,
												double* const state,
												const POPPARAM* const popparam,
												const double endtime,
												const double* rates)
{
	ADVANCER_THREECOMP* advanthreecomp = (ADVANCER_THREECOMP*)advan;	/* up cast */

	(void)record;
	(void)popparam;

	/* this advancer can only allow doses to first compartment */
	assert(advan->initcount > 0);
	assert(rates[0] >= 0.);
	assert(rates[1] == 0.);
	assert(rates[2] == 0.);

	/* reset the internal imodel IMODEL if necessary */
	if (advan->initcount != advanthreecomp->last_reset_count) {
		imodel_threecomp_reset_macro_constants(advanthreecomp, imodel);
		advanthreecomp->last_reset_count = advan->initcount;
	}

	const double time = advan->time;

	assert (endtime > time);

	/* ADVAN 3-compartment code obtained from: Abuhelwa AY, Foster DJ, Upton RN.
	 * ADVAN-style analytical solutions for common pharmacokinetic imodels.
	   Journal of pharmacological and toxicological methods. 2015 Jun 30;73:42-8. */
	const double k12 = advanthreecomp->k12;
	const double k21 = advanthreecomp->k21;
	const double k13 = advanthreecomp->k13;
	const double k31 = advanthreecomp->k31;
	const double E1 = advanthreecomp->E1;
	const double E2 = advanthreecomp->E2;
	const double E3 = advanthreecomp->E3;
	const double lambda1 = advanthreecomp->lambda1;
	const double lambda2 = advanthreecomp->lambda2;
	const double lambda3 = advanthreecomp->lambda3;

#pragma push_macro("A1")
#pragma push_macro("A2")
#pragma push_macro("A3")
#define A1 (state[0])
#define A2 (state[1])
#define A3 (state[2])
	const double t = endtime - time; /* change in time */
	const double A1last = A1;
	const double A2last = A2;
	const double A3last = A3;
	const double Doserate = rates[0];

	const double B = A2last*k21+A3last*k31;
	const double C = E3*A2last*k21+E2*A3last*k31;
	const double I = A1last*k12*E3-A2last*k13*k31+A3last*k12*k31;
	const double J = A1last*k13*E2+A2last*k13*k21-A3last*k12*k21;

	/* cached values for greater speed */
	const double exp_t_lambda1 = exp(-t*lambda1);
	const double exp_t_lambda2 = exp(-t*lambda2);
	const double exp_t_lambda3 = exp(-t*lambda3);

	const double A1term1 = A1last*(exp_t_lambda1*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp_t_lambda2*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp_t_lambda3*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
	const double A1term2 = exp_t_lambda1*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp_t_lambda2*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp_t_lambda3*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2));
	/* small optimization where large term does not have to be calulated in Doserate is zero */
	const double A1term3 = (Doserate == 0.) ? 0. : Doserate*((E2*E3)/(lambda1*lambda2*lambda3)-exp_t_lambda1*(E2-lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp_t_lambda2*(E2-lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp_t_lambda3*(E2-lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)));
	A1 = A1term1+A1term2+A1term3;    /* Amount in the central compartment */

	const double A2term1 = A2last*(exp_t_lambda1*(E1-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp_t_lambda2*(E1-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp_t_lambda3*(E1-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
	const double A2term2 = exp_t_lambda1*(I-A1last*k12*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp_t_lambda2*(A1last*k12*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp_t_lambda3*(A1last*k12*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2));
	/* small optimization where large term does not have to be calulated in Doserate is zero */
	const double A2term3 = (Doserate == 0.) ? 0. : Doserate*k12*(E3/(lambda1*lambda2*lambda3)-exp_t_lambda1*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp_t_lambda2*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp_t_lambda3*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)));
	A2 = A2term1+A2term2+A2term3;    /* Amount in the first-peripheral compartment */

	const double A3term1 = A3last*(exp_t_lambda1*(E1-lambda1)*(E2-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp_t_lambda2*(E1-lambda2)*(E2-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp_t_lambda3*(E1-lambda3)*(E2-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
	const double A3term2 = exp_t_lambda1*(J-A1last*k13*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp_t_lambda2*(A1last*k13*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp_t_lambda3*(A1last*k13*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2));
	/* small optimization where large term does not have to be calulated in Doserate is zero */
	const double A3term3 = (Doserate == 0.) ? 0. : Doserate*k13*(E2/(lambda1*lambda2*lambda3)-exp_t_lambda1*(E2-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp_t_lambda2*(E2-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp_t_lambda3*(E2-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)));
	A3 = A3term1+A3term2+A3term3;  /* Amount in the second-peripheral compartment */

#pragma pop_macro("A1")
#pragma pop_macro("A2")
#pragma pop_macro("A3")
}

ADVANFUNCS* pmx_advan_threecomp(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);
	assert(advanconfig->nstate == 0 || advanconfig->nstate == 3);

	let retinit = (ADVANTABLE_THREECOMP) {
		.advanfuncs = {
			.advan_size = sizeof(ADVANCER_THREECOMP),
			.construct = advancer_threecomp_construct,
			.destruct = advancer_threecomp_destruct,
			.info = advancer_threecomp_info,

			.reset = 0,
			.interval = advancer_threecomp_advance_interval,

			.advanconfig = advanconfig,
			.recordinfo = recordinfo_init(dataconfig),
			.nstate = 3,
		},
		.offsetV1 = structinfo_find_offset("V1", &advanconfig->imodelfields),
		.offsetV2 = structinfo_find_offset("V2", &advanconfig->imodelfields),
		.offsetV3 = structinfo_find_offset("V3", &advanconfig->imodelfields),
		.offsetCL = structinfo_find_offset("CL", &advanconfig->imodelfields),
		.offsetQ2 = structinfo_find_offset("Q2", &advanconfig->imodelfields),
		.offsetQ3 = structinfo_find_offset("Q3", &advanconfig->imodelfields),
	};
	assert(retinit.offsetV1 >= 0);
	assert(retinit.offsetV2 >= 0);
	assert(retinit.offsetV3 >= 0);
	assert(retinit.offsetCL >= 0);
	assert(retinit.offsetQ2 >= 0);
	assert(retinit.offsetQ3 >= 0);

	/* make binary copy so init can have const members */
	ADVANTABLE_THREECOMP* ret = malloc(sizeof(ADVANTABLE_THREECOMP));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANTABLE_THREECOMP));

	return &ret->advanfuncs;
}

