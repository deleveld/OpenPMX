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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "predict.h"
#include "ievaluate.h"
#include "idata.h"
#include "scatter.h"
#include "utils/c22.h"
#include "utils/various.h"
#include "pmxstate.h"

/* NOTE: this function must be thread safe on the level of an individual */
static void idata_predict_yhat_thread(INDIVID* const individ,
									  const ADVANFUNCS* const advanfuncs,
									  const POPMODEL* const popmodel,
									  const NONZERO* const nonzero,
									  const OPTIONS* const options,
									  const SCATTEROPTIONS* const scatteroptions)
{
	/* predict individual does not use reduced omega, this is just a sanity check */
	assert(nonzero == 0);
	(void)nonzero;
	(void)options;

	double etaarray[OPENPMX_OMEGA_MAX];
	forcount(i, OPENPMX_OMEGA_MAX)
		etaarray[i] = NAN;
	memcpy(etaarray, individ->eta, popmodel->nomega * sizeof(double));

	let ievaluate_args = ievaluate_args_init(individ->record,
											 individ->nrecord,
											 advanfuncs,
											 popmodel->theta,
											 popmodel->ntheta,
											 etaarray,
											 popmodel->nomega,
											 popmodel->sigma,
											 popmodel->nsigma,
											 scatteroptions ? scatteroptions->logstream : 0);
	struct timespec t3;
	clock_gettime(CLOCK_REALTIME, &t3);

	/* reset yhatvar when we write yhat */
	memset(individ->yhatvar, 0, individ->nrecord * sizeof(double));
	
	individual_evaluate(&ievaluate_args,
						individ->imodel,		/* do write imodel */
						individ->predictvars,	/* do write predictvars */
						individ->istate,		/* do write state */
						individ->yhat,			/* do write yhat */
						0, 						/* dont write yhatvar */
						0, 0);
	individ->ineval += 1;
	timespec_duration(&t3, &individ->eval_msec);
}

/* NOTE: this function must be thread safe on the level of an individual */
static void idata_predict_pred_thread(INDIVID* const individ,
									  const ADVANFUNCS* const advanfuncs,
									  const POPMODEL* const popmodel,
									  const NONZERO* const nonzero,
									  const OPTIONS* const options,
									  const SCATTEROPTIONS* const scatteroptions)
{
	/* predict individual does not use reduced omega, this is just a sanity check */
	assert(nonzero == 0);
	(void)nonzero;
	(void)options;

	/* zero eta for pred */
	double etaarray[OPENPMX_OMEGA_MAX];
	forcount(i, OPENPMX_OMEGA_MAX)
		etaarray[i] = NAN;
	forcount(i, popmodel->nomega)
		etaarray[i] = 0.;

	let ievaluate_args = ievaluate_args_init(individ->record,
											 individ->nrecord,
											 advanfuncs,
											 popmodel->theta,
											 popmodel->ntheta,
											 etaarray,
											 popmodel->nomega,
											 popmodel->sigma,
											 popmodel->nsigma,
											 scatteroptions ? scatteroptions->logstream : 0);
	struct timespec t3;
	clock_gettime(CLOCK_REALTIME, &t3);

	/* only write pred */
	individual_evaluate(&ievaluate_args,
						0,				/* dont write imodel */
						0,				/* dont write predictvars */
						0,				/* dont write state */
						individ->pred,	/* do write pred instead of yhat */
						0, 				/* dont write yhatvar */
						0, 0);
	individ->ineval += 1;
	timespec_duration(&t3, &individ->eval_msec);
}

static void idata_predict_yhat(const IDATA* const idata,
							   const ADVANFUNCS* const advanfuncs,
							   const POPMODEL* const popmodel,
							   const OPTIONS* const options)
{
	/* fill in yhat, imodel, and  predictvars */
	SCATTEROPTIONS scatteroptions = { };
	scatter_threads(idata, advanfuncs, popmodel, 0, options, &scatteroptions, idata_predict_yhat_thread);
}

void idata_predict_pred(const IDATA* const idata,
						const ADVANFUNCS* const advanfuncs,
						const POPMODEL* const popmodel,
						const OPTIONS* const options)
{
	/* fill in pred */
	SCATTEROPTIONS scatteroptions = { };
	scatter_threads(idata, advanfuncs, popmodel, 0, options, &scatteroptions, idata_predict_pred_thread);
}

void pmx_predict(OPENPMX* pmx)
{
	pmxstate_ensure(pmx);
	var pstate = pmx->state;

	var options = options_init(pmx);
	var popmodel = popmodel_init(pmx);

	idata_predict_yhat(&pstate->idata,
					  pstate->advanfuncs,
					  &popmodel,
					  &options);
}

void pmx_predict_pred(OPENPMX* pmx)
{
	pmxstate_ensure(pmx);
	var pstate = pmx->state;

	var options = options_init(pmx);
	var popmodel = popmodel_init(pmx);

	idata_predict_pred(&pstate->idata,
						  pstate->advanfuncs,
						  &popmodel,
						  &options);
}

