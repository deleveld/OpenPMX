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

/// This file implements an advancer for a three compartment mammilary 
/// model.
///
/// The equations were obtained from: Abuhelwa AY, Foster DJ, Upton RN.
/// ADVAN-style analytical solutions for common pharmacokinetic imodels.
/// Journal of pharmacological and toxicological methods. 2015 Jun 30;73:42-8.

#include <assert.h>
#include <math.h>

#include "advan/advan.h"
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
    var advanthreecomp = container_of(advan, ADVANCER_THREECOMP, advan);


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
	let imodeloffsets = container_of(advanthreecomp->advan.advanfuncs, ADVANTABLE_THREECOMP, advanfuncs);

	let V1 = *(const double*)(((char*)imodel) + imodeloffsets->offsetV1);
	let V2 = *(const double*)(((char*)imodel) + imodeloffsets->offsetV2);
	let V3 = *(const double*)(((char*)imodel) + imodeloffsets->offsetV3);
	let CL = *(const double*)(((char*)imodel) + imodeloffsets->offsetCL);
	let Q12 = *(const double*)(((char*)imodel) + imodeloffsets->offsetQ2);
	let Q13 = *(const double*)(((char*)imodel) + imodeloffsets->offsetQ3);
	let k10 = CL/V1;
	let k12 = Q12/V1;
	let k21 = k12*V1/V2;
	let k13 = Q13/V1;
	let k31 = k13*V1/V3;
	let k20 = 0.;
	let k30 = 0.;
	let E1 = k10+k12+k13;
	let E2 = k21+k20;
	let E3 = k31+k30;

	/* calculate hybrid rate constants */
	let a = E1+E2+E3;
	let b = E1*E2+E3*(E1+E2)-k12*k21-k13*k31;
	let c = E1*E2*E3-E3*k12*k21-E2*k13*k31;

	let m = (3.*b - (a*a))/3.;
	let n = (2.*(a*a*a) - 9.*a*b + 27.*c)/27.;
	let Q = ((n*n))/4. + ((m*m*m))/27.;

	let alpha = sqrt(-1.*Q);
	let beta = -1.*n/2.;
	let gamma = sqrt((beta*beta)+(alpha*alpha));
	let theta = atan2(alpha,beta);

	let pow_gamma_1_3 = cbrt(gamma);
	let lambda1 = a/3. + pow_gamma_1_3*(cos(theta/3.) + sqrt(3.)*sin(theta/3.));
	let lambda2 = a/3. + pow_gamma_1_3*(cos(theta/3.) - sqrt(3.)*sin(theta/3.));
	let lambda3 = a/3. -(2.*pow_gamma_1_3*cos(theta/3.));

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

__attribute__ ((hot))
static void advancer_threecomp_advance_interval(ADVAN* advan,
												const IMODEL* const imodel,
												const RECORD* const record,
												double* const state,
												const POPPARAM* const popparam,
												const double endtime,
												const double* rates)
{
    var advanthreecomp = container_of(advan, ADVANCER_THREECOMP, advan);

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

	let time = advan->time;

	assert (endtime > time);

	/* ADVAN 3-compartment code obtained from: Abuhelwa AY, Foster DJ, Upton RN.
	 * ADVAN-style analytical solutions for common pharmacokinetic imodels.
	   Journal of pharmacological and toxicological methods. 2015 Jun 30;73:42-8. */
	let k12 = advanthreecomp->k12;
	let k21 = advanthreecomp->k21;
	let k13 = advanthreecomp->k13;
	let k31 = advanthreecomp->k31;
	let E1 = advanthreecomp->E1;
	let E2 = advanthreecomp->E2;
	let E3 = advanthreecomp->E3;
	let lambda1 = advanthreecomp->lambda1;
	let lambda2 = advanthreecomp->lambda2;
	let lambda3 = advanthreecomp->lambda3;

#define A1 (state[0])
#define A2 (state[1])
#define A3 (state[2])
	let t = endtime - time; /* change in time */
	let A1last = A1;
	let A2last = A2;
	let A3last = A3;
	let Doserate = rates[0];

	let B = A2last*k21+A3last*k31;
	let C = E3*A2last*k21+E2*A3last*k31;
	let I = A1last*k12*E3-A2last*k13*k31+A3last*k12*k31;
	let J = A1last*k13*E2+A2last*k13*k21-A3last*k12*k21;

	/* cached values for greater speed */
	let exp_t_lambda1 = exp(-t*lambda1);
	let exp_t_lambda2 = exp(-t*lambda2);
	let exp_t_lambda3 = exp(-t*lambda3);
	let div_lambda2131 = 1./((lambda2-lambda1)*(lambda3-lambda1));
	let div_lambda1213 = 1./((lambda1-lambda2)*(lambda1-lambda3));
	let div_lambda1223 = 1./((lambda1-lambda2)*(lambda2-lambda3));
	let div_lambda1232 = 1./((lambda1-lambda2)*(lambda3-lambda2));
	let div_lambda1323 = 1./((lambda1-lambda3)*(lambda2-lambda3));
	let div_lambda1332 = 1./((lambda1-lambda3)*(lambda3-lambda2));
	let div_lambda12131 = 1./(lambda1*(lambda2-lambda1)*(lambda3-lambda1));
	let div_lambda21232 = 1./(lambda2*(lambda1-lambda2)*(lambda3-lambda2));
	let div_lambda31323 = 1./(lambda3*(lambda1-lambda3)*(lambda2-lambda3));
	let div_lambda123 = 1./(lambda1*lambda2*lambda3);
	let E1l1 = E1-lambda1;
	let E2l1 = E2-lambda1;
	let E3l1 = E3-lambda1;
	let E2l2 = E2-lambda2;
	let E2l3 = E2-lambda3;
	let E3l3 = E3-lambda3;
	let E3l2 = E3-lambda2;

	let A1term1 = A1last*(exp_t_lambda1*E2l1*E3l1*div_lambda2131
						  + exp_t_lambda2*E2l2*E3l2*div_lambda1232
						  + exp_t_lambda3*E2l3*E3l3*div_lambda1323);
	let A1term2 = exp_t_lambda1*(C-B*lambda1)*div_lambda1213
				+ exp_t_lambda2*(B*lambda2-C)*div_lambda1223
				+ exp_t_lambda3*(B*lambda3-C)*div_lambda1332;
	/* small optimization where large term does not have to be calulated in Doserate is zero */
	let A1term3 = (Doserate == 0.) ? 0. : Doserate*((E2*E3)*div_lambda123
													 - exp_t_lambda1*E2l1*E3l1*div_lambda12131
													 - exp_t_lambda2*E2l2*E3l2*div_lambda21232
													 - exp_t_lambda3*E2l3*E3l3*div_lambda31323);
	A1 = A1term1 + A1term2 + A1term3;    /* Amount in the central compartment */

	let A2term1 = A2last*(exp_t_lambda1*E1l1*E3l1*div_lambda2131
						  + exp_t_lambda2*(E1-lambda2)*E3l2*div_lambda1232
						  + exp_t_lambda3*(E1-lambda3)*E3l3*div_lambda1323);
	let A1last_k12 = A1last*k12;
	let A2term2 = exp_t_lambda1*(I-A1last_k12*lambda1)*div_lambda1213
				+ exp_t_lambda2*(A1last_k12*lambda2-I)*div_lambda1223
				+ exp_t_lambda3*(A1last_k12*lambda3-I)*div_lambda1332;
	/* small optimization where large term does not have to be calulated in Doserate is zero */
	let A2term3 = (Doserate == 0.) ? 0. : Doserate*k12*(E3*div_lambda123
														 - exp_t_lambda1*E3l1*div_lambda12131
														 - exp_t_lambda2*E3l2*div_lambda21232
														 - exp_t_lambda3*E3l3*div_lambda31323);
	A2 = A2term1 + A2term2 + A2term3;    /* Amount in the first-peripheral compartment */

	let A3term1 = A3last*(exp_t_lambda1*E1l1*E2l1*div_lambda2131
						  + exp_t_lambda2*(E1-lambda2)*E2l2*div_lambda1232
						  + exp_t_lambda3*(E1-lambda3)*E2l3*div_lambda1323);
	let A1last_k13 = A1last*k13;
	let A3term2 = exp_t_lambda1*(J-A1last_k13*lambda1)*div_lambda1213
				+ exp_t_lambda2*(A1last_k13*lambda2-J)*div_lambda1223
				+ exp_t_lambda3*(A1last_k13*lambda3-J)*div_lambda1332;
	/* small optimization where large term does not have to be calulated in Doserate is zero */
	let A3term3 = (Doserate == 0.) ? 0. : Doserate*k13*(E2*div_lambda123
														- exp_t_lambda1*E2l1*div_lambda12131
														- exp_t_lambda2*E2l2*div_lambda21232
														- exp_t_lambda3*E2l3*div_lambda31323);
	A3 = A3term1 + A3term2 + A3term3;  /* Amount in the second-peripheral compartment */
#undef A1
#undef A2
#undef A3
}

static inline void ensure(const int flag, const char* message)
{
	if (!flag) {
		printf("%s\n", message);
		exit(EXIT_FAILURE);
	}
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
	ensure(retinit.offsetV1 >= 0, "fatal: could not find V1");
	ensure(retinit.offsetV2 >= 0, "fatal: could not find V2");
	ensure(retinit.offsetV3 >= 0, "fatal: could not find V3");
	ensure(retinit.offsetCL >= 0, "fatal: could not find CL");
	ensure(retinit.offsetQ2 >= 0, "fatal: could not find Q2");
	ensure(retinit.offsetQ3 >= 0, "fatal: could not find Q3");

	/* make binary copy so init can have const members */
	ADVANTABLE_THREECOMP* ret = malloc(sizeof(ADVANTABLE_THREECOMP));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANTABLE_THREECOMP));

	return &ret->advanfuncs;
}

