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
	double V1, V2, V3, CL, Q2, Q3;

	/* parameter-only constants */
	double l[3];        /* lambda */
	double Mc[12];      /* time-independent part of M (infusion steady-state terms) */
	double K[3][12];    /* per-mode coefficients, already multiplied by 1/prod(l_i - l_j) */

	/* time-dependent cache: new_state[r] = M[r][0..2].state + M[r][3].rate */
	int tvalid;
	double tlast;
	double M[12];       /* 3x4 row-major */
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
} ADVANFUNCS_THREECOMP;

static void advancer_threecomp_setup(ADVANCER_THREECOMP* s,
									 const double V1, const double V2, const double V3,
									 const double CL, const double Q2, const double Q3)
{
	s->V1 = V1; s->V2 = V2; s->V3 = V3;
	s->CL = CL; s->Q2 = Q2; s->Q3 = Q3;

	/* division 1 of 2: 1/V1, 1/V2, 1/V3 */
	let iV = 1./(V1*V2*V3);
	let iV1 = V2*V3*iV;
	let iV2 = V1*V3*iV;
	let iV3 = V1*V2*iV;
	let k10 = CL*iV1;
	let k12 = Q2*iV1;
	let k21 = Q2*iV2;
	let k13 = Q3*iV1;
	let k31 = Q3*iV3;
	let E1 = k10 + k12 + k13;
	let E2 = k21;
	let E3 = k31;

	/* hybrid rate constants (trigonometric cubic solution) */
	let A = E1 + E2 + E3;
	let B = E1*E2 + E3*(E1+E2) - k12*k21 - k13*k31;
	let C = E1*E2*E3 - E3*k12*k21 - E2*k13*k31;

	let m = (3.*B - A*A)/3.;
	let n = (2.*A*A*A - 9.*A*B + 27.*C)/27.;
	let Q = n*n/4. + m*m*m/27.;

	let alpha = sqrt(-Q);
	let beta = -0.5*n;
	let gamma = sqrt(beta*beta + alpha*alpha);
	let theta = atan2(alpha, beta);

	let r = cbrt(gamma);
	double sn, cs;
	__builtin_sincos(theta/3., &sn, &cs);
	let s3 = 1.7320508075688772935274463415059*sn;

	let l1 = A/3. + r*(cs + s3);
	let l2 = A/3. + r*(cs - s3);
	let l3 = A/3. - 2.*r*cs;

	/* division 2 of 2: 1/(x12 x13 x23), 1/(l1 l2 l3), 1/l_i all from one reciprocal */
	let x12 = l1 - l2;
	let x13 = l1 - l3;
	let x23 = l2 - l3;
	let X = x12*x13*x23;
	let P = l1*l2*l3;
	let inv = 1./(X*P);
	let invden = inv*P;
	let ip = inv*X;			/* 1/(l1 l2 l3) */

	const double l[3] = { l1, l2, l3 };
	const double d[3] = { x23*invden, -x13*invden, x12*invden };
	const double il[3] = { l2*l3*ip, l1*l3*ip, l1*l2*ip };

	let k13k31 = k13*k31;
	let k12k31 = k12*k31;
	let k13k21 = k13*k21;
	let k12k21 = k12*k21;

	for (int i = 0; i < 3; i++) {
		let a = E1 - l[i];
		let b = E2 - l[i];
		let c = E3 - l[i];
		let di = d[i];
		let ili = il[i];
		double* const K = s->K[i];

		s->l[i] = l[i];

		K[0]  = di*b*c;               K[1]  = di*k21*c;            K[2]  = di*k31*b;            K[3]  = -di*ili*b*c;
		K[4]  = di*k12*c;             K[5]  = di*(c*a - k13k31);   K[6]  = di*k12k31;           K[7]  = -di*k12*ili*c;
		K[8]  = di*k13*b;             K[9]  = di*k13k21;           K[10] = di*(b*a - k12k21);   K[11] = -di*k13*ili*b;
	}

	for (int j = 0; j < 12; j++)
		s->Mc[j] = 0.;
	s->Mc[3]  = E2*E3*ip;
	s->Mc[7]  = k12*E3*ip;
	s->Mc[11] = k13*E2*ip;

	s->tvalid = 0;
	s->valid = 1;
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
	var s = container_of(advan, ADVANCER_THREECOMP, advan);

	(void)record;
	(void)popparam;

	assert(advan->initcount > 0);
	assert(endtime > advan->time);
	assert(rates[0] >= 0.);
	assert(rates[1] == 0.);
	assert(rates[2] == 0.);

	let o = container_of(advan->advanfuncs, ADVANFUNCS_THREECOMP, advanfuncs);
	let V1 = *(const double*)(((const char*)imodel) + o->offsetV1);
	let V2 = *(const double*)(((const char*)imodel) + o->offsetV2);
	let V3 = *(const double*)(((const char*)imodel) + o->offsetV3);
	let CL = *(const double*)(((const char*)imodel) + o->offsetCL);
	let Q2 = *(const double*)(((const char*)imodel) + o->offsetQ2);
	let Q3 = *(const double*)(((const char*)imodel) + o->offsetQ3);

	if (__builtin_expect(!s->valid || V1 != s->V1 || V2 != s->V2 || V3 != s->V3
						 || CL != s->CL || Q2 != s->Q2 || Q3 != s->Q3, 0))
		advancer_threecomp_setup(s, V1, V2, V3, CL, Q2, Q3);

	let t = endtime - advan->time;
	if (!s->tvalid || t != s->tlast) {
		let e0 = exp(-t*s->l[0]);
		let e1 = exp(-t*s->l[1]);
		let e2 = exp(-t*s->l[2]);
#pragma GCC unroll 12
		for (int j = 0; j < 12; j++)
			s->M[j] = s->Mc[j] + e0*s->K[0][j] + e1*s->K[1][j] + e2*s->K[2][j];
		s->tlast = t;
		s->tvalid = 1;
	}

	let A1 = state[0];
	let A2 = state[1];
	let A3 = state[2];
	let D  = rates[0];
	const double* const M = s->M;

	state[0] = M[0]*A1 + M[1]*A2 + M[2]*A3  + M[3]*D;
	state[1] = M[4]*A1 + M[5]*A2 + M[6]*A3  + M[7]*D;
	state[2] = M[8]*A1 + M[9]*A2 + M[10]*A3 + M[11]*D;
}

ADVANFUNCS* pmx_advan_threecomp(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);

	var nstate = advanconfig->nstate;
	if (nstate == 0)
		nstate = 3;
	int err = advan_ensure(nstate >= 3, __func__, "nstate should be >= 3");
	if (err)
		return 0;

	let retinit = (ADVANFUNCS_THREECOMP) {
		.advanfuncs = {
			.advan_size = sizeof(ADVANCER_THREECOMP),
			.construct = advancer_threecomp_construct,
			.destruct = advancer_threecomp_destruct,
			.info = advancer_threecomp_info,

			.reset = 0,
			.interval = advancer_threecomp_advance_interval,

			.advanconfig = advanconfig,
			.recordinfo = recordinfo_init(dataconfig),
			.nstate = nstate,
		},
		.offsetV1 = structinfo_find_offset("V1", &advanconfig->imodelfields),
		.offsetV2 = structinfo_find_offset("V2", &advanconfig->imodelfields),
		.offsetV3 = structinfo_find_offset("V3", &advanconfig->imodelfields),
		.offsetCL = structinfo_find_offset("CL", &advanconfig->imodelfields),
		.offsetQ2 = structinfo_find_offset("Q2", &advanconfig->imodelfields),
		.offsetQ3 = structinfo_find_offset("Q3", &advanconfig->imodelfields),
	};
	
	err |= advan_ensure(retinit.offsetV1 >= 0, __func__, "could not find V1");
	err |= advan_ensure(retinit.offsetV2 >= 0, __func__, "could not find V2");
	err |= advan_ensure(retinit.offsetV3 >= 0, __func__, "could not find V3");
	err |= advan_ensure(retinit.offsetCL >= 0, __func__, "could not find CL");
	err |= advan_ensure(retinit.offsetQ2 >= 0, __func__, "could not find Q2");
	err |= advan_ensure(retinit.offsetQ3 >= 0, __func__, "could not find Q3");
	if (err)
		return 0;

	/* make binary copy so init can have const members */
	ADVANFUNCS_THREECOMP* ret = malloc(sizeof(ADVANFUNCS_THREECOMP));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANFUNCS_THREECOMP));

	return &ret->advanfuncs;
}

