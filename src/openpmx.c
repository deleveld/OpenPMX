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
#include <values.h>
#include <unistd.h>

#include "openpmx.h"
#include "print.h"
#include "scatter.h"
#include "utils/c22.h"
#include "advan/advan.h"
#include "openpmx_internal.h"

#include "openpmx_compile_options.h"

STAGE1CONFIG stage1config_default(const STAGE1CONFIG* const stage1)
{
	STAGE1CONFIG ret = { 0 };
	if (stage1)
		ret = *stage1;

	let v = pow(DBL_EPSILON, 1./3.);
	if (ret.gradient_step == 0.)
		ret.gradient_step = 1e-4;
	if (ret.step_initial == 0.)
		ret.step_initial = 0.2;
	if (ret.step_refine == 0.)
		ret.step_refine = 0.1;
	if (ret.step_final == 0.)
		ret.step_final = v;
	if (ret.maxeval == 0)
		ret.maxeval = 1000;

	/* 	default omit_icov_resample = false */
	if (ret.icov_resample_tol == 0)
		ret.icov_resample_tol = DBL_MAX;

	return ret;
}

ESTIMCONFIG estimconfig_default(const ESTIMCONFIG* const estimate)
{
	ESTIMCONFIG ret = { 0 };
	if (estimate)
		ret = *estimate;

	ret.stage1 = stage1config_default(&ret.stage1);

	let v = pow(DBL_EPSILON, 1./3.);
	if (ret.optim.step_initial == 0.)
		ret.optim.step_initial = 0.2;
	if (ret.optim.step_refine == 0.)
		ret.optim.step_refine = 0.05;
	if (ret.optim.step_final == 0.)
		ret.optim.step_final = v * 10.;
	if (ret.optim.maxeval == 0)
		ret.optim.maxeval = 10000;

	if (ret.posthoc.deriv.step == 0.)
		ret.posthoc.deriv.step = ret.optim.step_initial;

	return ret;
}

SIMCONFIG simconfig_default(const SIMCONFIG* const simulate)
{
	SIMCONFIG ret = { 0 };
	if (simulate)
		ret = *simulate;

	if (ret.seed == 0)
		ret.seed = 200501041406; /* Birthdays of Joyce, Nette and Isabel :) */

	return ret;
}

OPTIONS options_default(const OPTIONS* const opt1)
{
	assert(opt1);
	OPTIONS ret = *opt1;

	/* nthreads == 0 means default number */
	if (ret.nthread == 0) {
		/* leave some cores over for the operating system */		
		int n = sysconf(_SC_NPROCESSORS_ONLN) - 4;
		if (n < 0)
			n = 0;
		ret.nthread = n;
	}

	ret.estimate = estimconfig_default(&ret.estimate);
	ret.simulate = simconfig_default(&ret.simulate);

	return ret;
}

POPPARAM popparam_init(const POPMODEL* const popmodel,
					   const ADVANFUNCS* const advanfuncs,
					   const double eta[static OPENPMX_OMEGA_MAX])
						{
	return (POPPARAM) {
		.theta = popmodel->theta,
		.ntheta = popmodel->ntheta,
		.eta = eta,
		.nomega = popmodel->nomega,
		.sigma = popmodel->sigma,
		.nsigma = popmodel->nsigma,
		.nstate = advanfuncs->nstate,
	};
}

OPTIONS options_init_from_pmx(const OPENPMX* const pmx)
{
	var ret = (OPTIONS) {
		.nthread = pmx->nthread,
		._offset1 = pmx->_offset1,

		.verbose = pmx->verbose,
		.brief = pmx->brief,
	};
	return options_default(&ret);
}

static PMXSTATE* pmxstate_alloc(const OPENPMX* const pmx)
{
	let advanfuncs = advanfuncs_alloc(&pmx->data, &pmx->advan);
	let popmodel = popmodel_init(pmx->theta, pmx->omega, pmx->sigma);

	var temp = (PMXSTATE) {
		.advanfuncs = advanfuncs,
		.idata = idata_construct(&advanfuncs->recordinfo,
								 popmodel.ntheta,
								 popmodel.nomega,
								 popmodel.nsigma,
								 advanfuncs->nstate,
								 advanfuncs->advanconfig->imodelfields.size,
								 advanfuncs->advanconfig->predictfields.size),
		.tablecount = 0,
		.rng = 0,
	};

	/* have to allocate because its a opaque pointer so we have to memcpy */
	var ret = mallocvar(PMXSTATE, 1);
	memcpy(ret, &temp, sizeof(temp));
	return ret;
}

void pmxstate_ensure(OPENPMX* const pmx)
{
	if (pmx->state == 0)
		pmx->state = pmxstate_alloc(pmx);
}

static void pmxstate_free(PMXSTATE* pstate)
{
	advanfuncs_free((void*)pstate->advanfuncs); /* keeping advanfuncs const in PMXSTATE is nice, but this cludge is needed to free */
	idata_destruct(&pstate->idata);

	if (pstate->rng)
		gsl_rng_free(pstate->rng);

	scatter_cleanup();
	free(pstate);
}

OPENPMX pmx_copy(const OPENPMX* const pmx)
{
	OPENPMX ret = *pmx;
	ret.state = 0;
	return ret;
}

void pmx_copy_popparams(OPENPMX* dest, const OPENPMX* const src)
{
	memcpy(dest->theta, src->theta, sizeof(dest->theta));
	memcpy(dest->sigma, src->sigma, sizeof(dest->sigma));
	memcpy(dest->omega, src->omega, sizeof(dest->omega));

	dest->result = (PMXRESULT) { .objfn = DBL_MAX,
								 .type = OBJFN_INVALID,
								 .nfunc = 0 };
}

void pmx_cleanup(OPENPMX* pmx)
{
	pmxstate_free(pmx->state);
	pmx->state = 0;
}


