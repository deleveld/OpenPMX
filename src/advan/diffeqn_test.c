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
	ADVANCER_DIFFEQN_CALLBACK_ARGS args;
	double hstep;
} ADVANCER_TESTRK4;

typedef struct {
	ADVANFUNCS advanfuncs;
	ADVAN_DIFFEQN diffeqn;
} ADVANTABLE_TESTRK4;

static inline double get_hstep(const struct ADVANCONFIG* const advanconfig)
{
	var hstep = advanconfig->args.diffeqn.hstart;
	if (hstep == 0.)
		hstep = 10./60.;
	return hstep;
}

static void advancer_diffeqn_info(const struct ADVANFUNCS* const advanfuncs, FILE* f)
{
	fprintf(f, "advan model: DES RK4 test solver\n");
	fprintf(f, "advan nstate: %i\n", advanfuncs->nstate);
	fprintf(f, "advan hstep: %g\n", get_hstep(advanfuncs->advanconfig));
}

__attribute__ ((hot))
static void advancer_diffeqn_wrapper(double TIME, const double* const A, double* DADT, const ADVANCER_DIFFEQN_CALLBACK_ARGS *args)
{
	let diffeqn = args->diffeqn;
	let imodel = args->imodel;
	let record = args->record;
	let popparam = args->popparam;
	let rates = args->rates;
	let nstate = args->nstate;

	/* pass imodel down to user DES */
	diffeqn(DADT, imodel, record, A, popparam, TIME);

	/* after user DES we add the infusion rates, so the user does not see these */
	assert(rates);
	forcount(i, nstate)
		DADT[i] += rates[i];
}

static void advancer_diffeqn_construct(ADVAN* advan, const ADVANFUNCS* const advanfuncs)
{
	var advandes = (ADVANCER_TESTRK4*)advan; /* cast up */
	let destest = (ADVANTABLE_TESTRK4*)advanfuncs;

	assert(advanfuncs->advan_size == sizeof(ADVANCER_TESTRK4));
	advan_base_construct(advan, advanfuncs); /* zeros full size */

	/* we only need to set these once */
	/* the rest need to be set at each interval */
	advandes->args.diffeqn = destest->diffeqn;
	advandes->args.nstate = advanfuncs->nstate;

	advandes->hstep = get_hstep(advanfuncs->advanconfig);
}

static void advancer_diffeqn_destruct(ADVAN * advan)
{
	advan_base_destruct(advan);
}

static void advancer_diffeqn_reset(ADVAN * advan, const int full)
{
	(void) advan;
	(void) full;
}

static void advancer_diffeqn_advance_interval(ADVAN* advan,
											  const IMODEL* const imodel,
											  const RECORD* const record,
											  double* const state,
											  const POPPARAM* const popparam,
											  const double endtime,
											  const double* rates)
{
	ADVANCER_TESTRK4* advandes = (ADVANCER_TESTRK4*)advan;	/* up cast */
	/* advandes->args.diffeqn = already set in constructor */
	advandes->args.advan = advan;
	advandes->args.imodel = imodel;
	advandes->args.record = record;
	advandes->args.popparam = popparam;
	advandes->args.rates = rates;
	/* advandes->args.nstate = nstate; already set in constructor */

	assert(advan->initcount > 0);

	let nstate = advandes->args.nstate;
	double k1[nstate];
	double k2[nstate];
	double k3[nstate];
	double k4[nstate];
	double teststate[nstate];

	while (advan->time < endtime) {

		/* upper limit to the step size we can take */
		var stepend = advan->time + advandes->hstep;
		if (stepend > endtime)
			stepend = endtime;
		let stepsize = stepend - advan->time;

		advancer_diffeqn_wrapper(advan->time, state, k1, &advandes->args);

		forcount(i, nstate)
			teststate[i] = state[i] + stepsize * k1[i] / 2.;
		advancer_diffeqn_wrapper(advan->time + stepsize / 2., teststate, k2, &advandes->args);

		forcount(i, nstate)
			teststate[i] = state[i] + stepsize * k2[i] / 2.;
		advancer_diffeqn_wrapper(advan->time + stepsize / 2., teststate, k3, &advandes->args);

		forcount(i, nstate)
			teststate[i] = state[i] + stepsize * k3[i];
		advancer_diffeqn_wrapper(advan->time + stepsize, teststate, k4, &advandes->args);

		forcount(i, nstate)
			state[i] += stepsize * (k1[i] + 2.*k2[i] + 2.*k3[i] + k4[i]) / 6.;

		/* go to the end of the interval */
		advan->time = stepend;

	}
	assert(advan->time == endtime); /* is this redundant?, the ode solver may have already set this */
}

ADVANFUNCS* pmx_advan_diffeqn_test(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);
	assert(advanconfig->diffeqn);
	assert(advanconfig->nstate);

	let retinit = (ADVANTABLE_TESTRK4) {
		.advanfuncs = {
			.advan_size = sizeof(ADVANCER_TESTRK4),
			.construct = advancer_diffeqn_construct,
			.destruct = advancer_diffeqn_destruct,
			.info = advancer_diffeqn_info,

			.reset = advancer_diffeqn_reset,
			.interval = advancer_diffeqn_advance_interval,

			.advanconfig = advanconfig,
			.recordinfo = recordinfo_init(dataconfig),
			.nstate = advanconfig->nstate,
		},
		.diffeqn = advanconfig->diffeqn,
	};

	ADVANTABLE_TESTRK4* ret = malloc(sizeof(ADVANTABLE_TESTRK4));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANTABLE_TESTRK4));

	return &ret->advanfuncs;
}

