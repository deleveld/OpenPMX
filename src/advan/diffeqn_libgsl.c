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

/// This file implements the ODE solvers inplemented in libGSL.

#include <assert.h>
#include <math.h>
#include <string.h>

#include "advan/advan.h"
#include "utils/c22.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

typedef struct {
	ADVAN advan;

	gsl_odeiv2_system sys;
	gsl_odeiv2_driver* d;

	ADVANCER_DIFFEQN_CALLBACK_ARGS args;
} ADVANCER_LIBGSL;

#define STEPTYPE_NAME_LENGTH	64
#define DEFAULT_ABSTOL			1e-9
#define DEFAULT_RELTOL			1e-9

/* TODO: Should these be moved to ADVANCER_LIBGSL so that changes to the config actually get passed on when an ADVAN is made */
typedef struct {
	ADVANFUNCS advanfuncs;
	ADVAN_DIFFEQN diffeqn;
	double abstol;
	double reltol;
	double hstart;
	const gsl_odeiv2_step_type* steptype;
	char steptype_name[STEPTYPE_NAME_LENGTH];
} ADVANTABLE_LIBGSL;

static void advancer_diffeqn_libgsl_info(const struct ADVANFUNCS* const advanfuncs, FILE* f)
{
	let libgsl = (const ADVANTABLE_LIBGSL*)advanfuncs;

	fprintf(f, "advan model GNU Scientific Library ODE\n");
	fprintf(f, "advan stepper %s\n", libgsl->steptype_name);
	fprintf(f, "advan nstate %i\n", advanfuncs->nstate);
	fprintf(f, "advan abstol %g reltol %g hstart %g\n", libgsl->abstol, libgsl->reltol, libgsl->hstart);
}

__attribute__ ((hot))
static int advancer_diffeqn_libgsl_wrapper(double T, const double* A, double* DADT, void *params)
{
	let args = (const ADVANCER_DIFFEQN_CALLBACK_ARGS*)params;
	let diffeqn = args->diffeqn;
	let imodel = args->imodel;
	let record = args->record;
	let popparam = args->popparam;
	let rates = args->rates;
	let nstate = args->nstate;

	/* pass imodel down to user DES */
	memset(DADT, 0, nstate * sizeof(double));
	diffeqn(DADT, imodel, record, A, popparam, T);

	/* after user DES we add the infusion rates, so the user does not see these */
	assert(rates);
	forcount(i, nstate)
		DADT[i] += rates[i];

	return GSL_SUCCESS;
}

static void advancer_diffeqn_libgsl_construct(ADVAN* advan, const ADVANFUNCS* const advanfuncs)
{
	var advandes = (ADVANCER_LIBGSL*)advan; /* cast up */
	let libgsl = (ADVANTABLE_LIBGSL*)advanfuncs;

	assert(advanfuncs->advan_size == sizeof(ADVANCER_LIBGSL));
	advan_base_construct(advan, advanfuncs); /* zeros full size */

	let nstate = advanfuncs->nstate;
	advandes->sys = (gsl_odeiv2_system) {
		.function = advancer_diffeqn_libgsl_wrapper,
		.jacobian = NULL,
		.dimension = nstate,
		.params = &advandes->args,
	};

    advandes->d = gsl_odeiv2_driver_alloc_y_new(&advandes->sys, libgsl->steptype, libgsl->hstart, libgsl->abstol, libgsl->reltol);

	/* we only need to set these once */
	/* the rest need to be set at each interval */
	advandes->args.diffeqn = libgsl->diffeqn;
	advandes->args.nstate = nstate;
}

static void advancer_diffeqn_libgsl_destruct(ADVAN * advan)
{
	var advandes = (ADVANCER_LIBGSL*)advan; /* cast up */

	gsl_odeiv2_driver_free(advandes->d);

	advan_base_destruct(advan);
}

static void advancer_diffeqn_libgsl_reset(ADVAN * advan, const int full)
{
	var advandes = (ADVANCER_LIBGSL*)advan; /* cast up */
	let libgsl = (ADVANTABLE_LIBGSL*)advan->advanfuncs;

	gsl_odeiv2_driver_reset(advandes->d);
	if (full)
		gsl_odeiv2_driver_reset_hstart(advandes->d, libgsl->hstart);
}

static void advancer_diffeqn_libgsl_advance_interval(ADVAN* advan,
													 const IMODEL* const imodel,
													 const RECORD* const record,
													 double* const state,
													 const POPPARAM* const popparam,
													 const double endtime,
													 const double* rates)
{
	var advandes = (ADVANCER_LIBGSL*)advan;	/* up cast */

	/* advandes->args.diffeqn already set in constructor */
	advandes->args.advan = advan;
	advandes->args.imodel = imodel;
	advandes->args.record = record;
	advandes->args.popparam = popparam;
	advandes->args.rates = rates;
	/* advandes->args.nstate = nstate; already set in constructor */

	assert(advan->initcount > 0);

	while (advan->time < endtime) {
		let status = gsl_odeiv2_driver_apply(advandes->d, &advan->time, endtime, state);
		if (status != GSL_SUCCESS) {
			fprintf(stderr, "error integrating from %f to %f\n", advan->time, endtime);
			fprintf(stderr, " - delta is %f\n", endtime - advan->time);
			var status_string = "unknown";
			switch (status) {
				case GSL_FAILURE: status_string = "GSL_FAILURE"; break;
				case GSL_EMAXITER: status_string = "GSL_EMAXITER"; break;
				case GSL_ENOPROG: status_string = "GSL_ENOPROG"; break;
				case GSL_EBADFUNC: status_string = "GSL_EBADFUNC"; break;
				case GSL_SUCCESS: status_string = "GSL_SUCCESS"; break;
				case GSL_EFAULT: status_string = "GSL_EFAULT"; break;
				case GSL_EINVAL: status_string = "GSL_EINVAL"; break;
				default: break;
			}
			fprintf(stderr, " - status %s\n", status_string);
			exit(EXIT_FAILURE);
		}
	}
	assert(advan->time == endtime); /* is this redundant?, the ode solver may have already set this */
}

ADVANFUNCS* pmx_advan_diffeqn_libgsl(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->init);
	assert(advanconfig->predict);
	assert(advanconfig->diffeqn);
	assert(advanconfig->nstate);

	var retinit = (ADVANTABLE_LIBGSL) {
		.advanfuncs = {
			.advan_size = sizeof(ADVANCER_LIBGSL),
			.construct = advancer_diffeqn_libgsl_construct,
			.destruct = advancer_diffeqn_libgsl_destruct,
			.info = advancer_diffeqn_libgsl_info,

			.reset = advancer_diffeqn_libgsl_reset,
			.interval = advancer_diffeqn_libgsl_advance_interval,

			.advanconfig = advanconfig,
			.recordinfo = recordinfo_init(dataconfig),
			.nstate = advanconfig->nstate,
		},
		.diffeqn = advanconfig->diffeqn,
	};
	retinit.hstart = advanconfig->args.diffeqn.hstart;
	retinit.abstol = advanconfig->args.diffeqn.abstol;
	retinit.reltol = advanconfig->args.diffeqn.reltol;
	memset(retinit.steptype_name, 0, sizeof(retinit.steptype_name));

/// The default setting are:
/// - ODE stepping function `gsl_odeiv2_step_rk8pd`
	let steptype = advanconfig->args.diffeqn.steptype;
	var steptype_name = "";
	
	if (steptype == 0 || strcmp(steptype, "rk8pd") == 0) { 
		retinit.steptype = gsl_odeiv2_step_rk8pd;
		steptype_name = "rk8pd, Runge-Kutta Prince-Dormand (8,9)";

	} else if (strcmp(steptype, "msadams") == 0) {
		retinit.steptype = gsl_odeiv2_step_msadams;
		steptype_name = "msadams, Multistep Adams";

	} else if (strcmp(steptype, "rkf45") == 0) {
		retinit.steptype = gsl_odeiv2_step_rkf45;
		steptype_name = "rk45, Runge-Kutta-Fehlberg (4,5)";

	} else if (strcmp(steptype, "rkck") == 0) {
		retinit.steptype = gsl_odeiv2_step_rkck;
		steptype_name = "rkck, Runge-Kutta Cash-Karp (4,5) ";

	} else if (strcmp(steptype, "rk4") == 0) {
		retinit.steptype = gsl_odeiv2_step_rk4;
		steptype_name = "rk4, Explicit 4th order (classical) Runge-Kutta";

	} else if (strcmp(steptype, "rk2") == 0) {
		retinit.steptype = gsl_odeiv2_step_rk2;
		steptype_name = "rk2, Explicit embedded Runge-Kutta (2,3)";

	} else {
		fprintf(stderr, "error: Invalid steptype \"%s\"\n", steptype);
		exit(EXIT_FAILURE);
	}
	snprintf(retinit.steptype_name, STEPTYPE_NAME_LENGTH, "%s", steptype_name);

	if (retinit.hstart == 0.)
		retinit.hstart = 1.;
	if (retinit.abstol == 0.)
		retinit.abstol = DEFAULT_ABSTOL;
	if (retinit.reltol == 0.)
		retinit.reltol = DEFAULT_RELTOL;

	ADVANTABLE_LIBGSL* ret = malloc(sizeof(ADVANTABLE_LIBGSL));
	assert(ret);
	memcpy(ret, &retinit, sizeof(ADVANTABLE_LIBGSL));

	return &ret->advanfuncs;
}

