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
 
/// This file implements an advancer for a two compartment mammilary 
/// model.
///
/// The equations were obtained from: Abuhelwa AY, Foster DJ, Upton RN.
/// ADVAN-style analytical solutions for common pharmacokinetic imodels.
/// Journal of pharmacological and toxicological methods. 2015 Jun 30;73:42-8.
///
/// Applied some AI guided optimizations for speed.

#include <assert.h>
#include <math.h>
#include <string.h>

#include "advan/advan.h"
#include "utils/c22.h"

typedef struct {
	ADVAN advan;

	/* cache key: the parameters the constants below were built from */
	int valid;
	double V1, V2, CL, Q2;

	/* parameter-only constants */
	double l[2];        /* lambda */
	double Mc[6];       /* time-independent part of M */
	double K[2][6];     /* per-mode coefficients, already multiplied by 1/(l1-l2) */

	/* time-dependent cache: new_state[r] = M[r][0..1].state + M[r][2].rate */
	int tvalid;
	double tlast;
	double M[6];        /* 2x3 row-major */
} ADVANCER_TWOCOMP;

static void advancer_twocomp_info(const struct ADVANFUNCS* const advanfuncs, FILE* f)
{
	(void) advanfuncs;

	fprintf(f, "advan model two compartment\n");
	fprintf(f, "advan parameters V1 V2 CL Q2\n");
}
static void advancer_twocomp_construct(ADVAN* advan, const ADVANFUNCS* const advanfuncs)
{
    var advantwocomp = container_of(advan, ADVANCER_TWOCOMP, advan);

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
} ADVANFUNCS_TWOCOMP;

static void advancer_twocomp_setup(ADVANCER_TWOCOMP* s,
								   const double V1, const double V2,
								   const double CL, const double Q12)
{
	s->V1 = V1; s->V2 = V2; s->CL = CL; s->Q2 = Q12;

	/* division 1 of 2: 1/V1 and 1/V2 */
	let iV = 1./(V1*V2);
	let iV1 = V2*iV;
	let iV2 = V1*iV;
	let k12 = Q12*iV1;
	let k21 = Q12*iV2;
	let E1 = CL*iV1 + k12;
	let E2 = k21;

	/* hybrid rate constants */
	let sum = E1 + E2;
	let sd = sqrt(sum*sum - 4.*(E1*E2 - k12*k21));	/* lambda1 - lambda2 */
	let l1 = 0.5*(sum + sd);
	let l2 = 0.5*(sum - sd);

	/* division 2 of 2: 1/(l1-l2) and 1/(l1 l2) from one reciprocal */
	let l12 = l1*l2;
	let inv = 1./(sd*l12);
	let dl = inv*l12;		/* 1/(l1-l2) */
	let p = inv*sd;			/* 1/(l1 l2) */
	let il1 = l2*p;
	let il2 = l1*p;

	s->l[0] = l1;
	s->l[1] = l2;

	double* const K1 = s->K[0];
	K1[0] = -dl*(E2 - l1);	K1[1] = -dl*k21;		K1[2] = -dl*(1. - E2*il1);
	K1[3] = -dl*k12;		K1[4] = -dl*(E1 - l1);	K1[5] =  dl*k12*il1;

	double* const K2 = s->K[1];
	K2[0] =  dl*(E2 - l2);	K2[1] =  dl*k21;		K2[2] =  dl*(1. - E2*il2);
	K2[3] =  dl*k12;		K2[4] =  dl*(E1 - l2);	K2[5] = -dl*k12*il2;

	s->Mc[0] = 0.; s->Mc[1] = 0.; s->Mc[2] = E2*p;
	s->Mc[3] = 0.; s->Mc[4] = 0.; s->Mc[5] = k12*p;

	s->tvalid = 0;
	s->valid = 1;
}

__attribute__ ((hot))
static void advancer_twocomp_advance_interval(ADVAN* advan,
											  const IMODEL* const imodel,
											  const RECORD* const record,
											  double* const state,
											  const POPPARAM* const popparam,
											  const double endtime,
											  const double* rates)
{
	var s = container_of(advan, ADVANCER_TWOCOMP, advan);

	(void)record;
	(void)popparam;

	assert(advan->initcount > 0);

	/* this advancer can only allow doses to first compartment */
	assert(endtime > advan->time);
	assert(rates[0] >= 0.);
	assert(rates[1] == 0.);

	let o = container_of(advan->advanfuncs, ADVANFUNCS_TWOCOMP, advanfuncs);
	let V1 = *(const double*)(((const char*)imodel) + o->offsetV1);
	let V2 = *(const double*)(((const char*)imodel) + o->offsetV2);
	let CL = *(const double*)(((const char*)imodel) + o->offsetCL);
	let Q12 = *(const double*)(((const char*)imodel) + o->offsetQ2);

	if (__builtin_expect(!s->valid || V1 != s->V1 || V2 != s->V2
						 || CL != s->CL || Q12 != s->Q2, 0))
		advancer_twocomp_setup(s, V1, V2, CL, Q12);

	let t = endtime - advan->time;
	if (!s->tvalid || t != s->tlast) {
		let e0 = exp(-t*s->l[0]);
		let e1 = exp(-t*s->l[1]);
#pragma GCC unroll 6
		for (int j = 0; j < 6; j++)
			s->M[j] = s->Mc[j] + e0*s->K[0][j] + e1*s->K[1][j];
		s->tlast = t;
		s->tvalid = 1;
	}

	let A1 = state[0];
	let A2 = state[1];
	let D  = rates[0];
	const double* const M = s->M;

	state[0] = M[0]*A1 + M[1]*A2 + M[2]*D;
	state[1] = M[3]*A1 + M[4]*A2 + M[5]*D;
}

ADVANFUNCS* pmx_advan_twocomp(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);

	var nstate = advanconfig->nstate;
	if (nstate == 0)
		nstate = 2;
	int err = advan_ensure(nstate >= 2, __func__, "nstate should be >= 2");
	if (err)
		return 0;

	let retinit = (ADVANFUNCS_TWOCOMP) {
		.advanfuncs = {
			.advan_size = sizeof(ADVANCER_TWOCOMP),
			.construct = advancer_twocomp_construct,
			.destruct = advancer_twocomp_destruct,
			.info = advancer_twocomp_info,

			.reset = 0,
			.interval = advancer_twocomp_advance_interval,

			.advanconfig = advanconfig,
			.recordinfo = recordinfo_init(dataconfig),
			.nstate = nstate,
		},
		.offsetV1 = structinfo_find_offset("V1", &advanconfig->imodelfields),
		.offsetV2 = structinfo_find_offset("V2", &advanconfig->imodelfields),
		.offsetCL = structinfo_find_offset("CL", &advanconfig->imodelfields),
		.offsetQ2 = structinfo_find_offset("Q2", &advanconfig->imodelfields),
	};
	
	err |= advan_ensure(retinit.offsetV1 >= 0, __func__, "could not find V1");
	err |= advan_ensure(retinit.offsetV2 >= 0, __func__, "could not find V2");
	err |= advan_ensure(retinit.offsetCL >= 0, __func__, "could not find CL");
	err |= advan_ensure(retinit.offsetQ2 >= 0, __func__, "could not find Q2");
	if (err)
		return 0;

	ADVANFUNCS_TWOCOMP* ret = malloc(sizeof(ADVANFUNCS_TWOCOMP));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANFUNCS_TWOCOMP));

	return &ret->advanfuncs;
}


