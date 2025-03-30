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
#include <float.h>

#include "ievaluate.h"
#include "print.h"
#include "advan/advan.h"
#include "dataconfig/dataconfig.h"
#include "dataconfig/recordinfo.h"
#include "utils/c22.h"

#include "openpmx_compile_options.h"

/* NOTE: these functions must be thread safe on the level of an individual */

__attribute__ ((hot))
static inline double evaluate_yhat(const IMODEL* const imodel,
								   const PREDICTSTATE* const predictstate,
								   const POPPARAM* const popparam,
								   const IMODEL_PREDICT predict,
								   const double errarray[static OPENPMX_SIGMA_MAX],
								   PREDICTVARS* predictvars)
{
	memset(predictvars, 0, OPENPMX_PREDICTVARS_MAX * sizeof(double));
		
	/* the errarray should be already set to zero for a normal call but
	 * it could be non-zero if we are simulating with residual error */
	return predict(imodel, predictstate, popparam, errarray, predictvars);
}

__attribute__ ((hot))
static double evaluate_yhatvar(const IMODEL* const imodel,
							   const PREDICTSTATE* const predictstate,
							   const POPPARAM* const popparam,
							   const IMODEL_PREDICT predict,
							   double errarray[static OPENPMX_SIGMA_MAX],
							   PREDICTVARS* predictvars)
{
	/* the errarray should be already set to zero and we preserve this
	 * across calls to this function */

	/* Do error propagation to get the variance of the prediction
	 * https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Simplification */
	var yhatvar = 0.;
	let sigma = popparam->sigma;
	let nsigma = popparam->nsigma;
	forcount(j, nsigma) {
		if (sigma[j] != 0.) {
			let g = sqrt(sigma[j]);

			let above = g;
			errarray[j] = above;
			let ya1 = evaluate_yhat(imodel, predictstate, popparam, predict, errarray, predictvars);

			let below = -g;
			errarray[j] = below;
			let ya2 = evaluate_yhat(imodel, predictstate, popparam, predict, errarray, predictvars);

			errarray[j] = 0;
			var deriv = (ya1 - ya2) / (above - below);

			yhatvar += (deriv*deriv) * sigma[j];
		}
	}
	return yhatvar;
}

static int assure_state_finite(const double* a, const int n)
{
	forcount(i, n) {
		if (isfinite(a[i]) != 1)
			return 0;
	}
	return 1;
}

typedef struct {
	double errarray[OPENPMX_SIGMA_MAX];
	double _predictvars[OPENPMX_PREDICTVARS_MAX];	/* will be cast to PREDICTVARS */
	double _imodel[OPENPMX_IMODEL_MAX];				/* will be cast to IMODEL */
} ADVAN_MODEL_MEMORY;

/* This is the core function evaluating an individual predictions to calculate
 * the objective function. Speeding up this function is quite important. */
 __attribute__ ((hot))
void individual_fasteval(const IEVALUATE_ARGS* const ievaluate_args,
						 double* const obs_lndet,
						 double* const obs_min2ll)
{
	let advanfuncs = ievaluate_args->advanfuncs;
	let record = ievaluate_args->record;
	let nrecord = ievaluate_args->nrecord;
	let popparam = &ievaluate_args->popparam;

	char advan_memory[advanfuncs->advan_size];
	var advan = (ADVAN*)advan_memory;
	advanfuncs->construct(advan, advanfuncs);

	ADVAN_MODEL_MEMORY advanmem = { 0 };
	var predictvars = (PREDICTVARS*)advanmem._predictvars;
	var imodel = (IMODEL*)advanmem._imodel;

	/* iterate over individuals data */
	let advanconfig = advanfuncs->advanconfig;
	let predict = advanconfig->predict;
	let recordinfo = &advanfuncs->recordinfo;
	let recordsize = recordinfo->dataconfig->recordfields.size;
	var local_obs_min2ll = 0.;
	var local_obs_lndet = 0.;
	const RECORD* ptr = record;
	forcount(i, nrecord) {
		let predictstate = advan_advance(advan, imodel, ptr, popparam);
		let evid = RECORDINFO_EVID(recordinfo, ptr);
		if (evid == 0) {
			let dv = RECORDINFO_DV(recordinfo, ptr);
			let yhat = evaluate_yhat(imodel, &predictstate, popparam, predict, advanmem.errarray, predictvars);
			let yhatvar = evaluate_yhatvar(imodel, &predictstate, popparam, predict, advanmem.errarray, predictvars);

			let err = dv - yhat;
			local_obs_min2ll += (err * err) / yhatvar;

			local_obs_lndet += log(yhatvar);
		}
		ptr = RECORD_INDEX(ptr, recordsize, 1);	
	}
	*obs_min2ll = local_obs_min2ll;
	*obs_lndet = local_obs_lndet;
		
	advanfuncs->destruct(advan);
}

/* This is a core function evaluating an individual by advancing over
 * the records and evaluating the imodel, predictions and objective function.
 * Speeding up this function is quite important.
 * I dont think making arguments restrict would help because the only pointer
 * read from has a different type than all of the other pointers, and these
 * are all only written to. */
 __attribute__ ((hot))
void individual_evaluate(const IEVALUATE_ARGS* const ievaluate_args,
						 IMODEL* const imodel_saved,
						 PREDICTVARS* const predictvars_saved,
						 double* const istate,
						 double* const YHAT,
						 double* const YHATVAR,
						 double* const obs_lndet,
						 double* const obs_min2ll)
{
	let advanfuncs = ievaluate_args->advanfuncs;
	let record = ievaluate_args->record;
	let nrecord = ievaluate_args->nrecord;
	let popparam = &ievaluate_args->popparam;

	char advan_memory[advanfuncs->advan_size];
	var advan = (ADVAN*)advan_memory;
	advanfuncs->construct(advan, advanfuncs);

	ADVAN_MODEL_MEMORY advanmem = { 0 };
	var predictvars = (PREDICTVARS*)advanmem._predictvars;
	var imodel = (IMODEL*)advanmem._imodel;

	/* iterate over individuals data */
	let advanconfig = advanfuncs->advanconfig;
	let predict = advanconfig->predict;
	let predictall = advanconfig->predictall;
	let nstate = advanfuncs->nstate;
	let recordinfo = &advanfuncs->recordinfo;
	let imodel_size = advanconfig->imodelfields.size;
	let predictvars_size = advanconfig->predictfields.size;
	var local_obs_min2ll = 0.;
	var local_obs_lndet = 0.;
	const RECORD* ptr = record;
	forcount(i, nrecord) {
		let predictstate = advan_advance(advan, imodel, ptr, popparam);

		let evid = RECORDINFO_EVID(recordinfo, ptr);
		var yhat = 0.;
		if (evid == 0 || predictall) 
			yhat = evaluate_yhat(imodel, &predictstate, popparam, predict, advanmem.errarray, predictvars);
		if (YHAT)
			YHAT[i] = yhat;

		/* save imodel, predictvars and state */
		if (imodel_saved) {
			let imodelptr = (IMODEL*)((char*)imodel_saved + i * imodel_size);
			memcpy(imodelptr, imodel, imodel_size);

			let predictvarsptr = (PREDICTVARS*)((char*)predictvars_saved + i * predictvars_size);
			memcpy(predictvarsptr, predictvars, predictvars_size);

			memcpy(&istate[i * nstate], predictstate.state, nstate * sizeof(double));
		}
		
		/* objective function only for observations */
		if (evid == 0) {
			let yhatvar = evaluate_yhatvar(imodel, &predictstate, popparam, predict, advanmem.errarray, predictvars);
			if (YHATVAR)
				YHATVAR[i] = yhatvar;

			/* we dont really need to have the if here but its probably faster if its in */
			let dv = RECORDINFO_DV(recordinfo, ptr);
			if (obs_min2ll) {
				let err = dv - yhat;
				let min2ll_term = (err * err) / yhatvar;
				local_obs_min2ll += min2ll_term;
			}
			if (obs_lndet) {
				let lndet_term = log(yhatvar);
				local_obs_lndet += lndet_term;
			}
		} else {
			/* yhatvar must be set to zero for non-observations since these will be used to calculate the
			 * individual covariance matrix */
			if (YHATVAR)
				YHATVAR[i] = 0.;
		}
		ptr = RECORDINFO_INDEX(recordinfo, ptr, 1);
	}
	if (obs_min2ll)
		*obs_min2ll = local_obs_min2ll;
	if (obs_lndet) 
		*obs_lndet = local_obs_lndet;
		
	advanfuncs->destruct(advan);
}

void individual_checkout(const IEVALUATE_ARGS* const ievaluate_args)
{
	let advanfuncs = ievaluate_args->advanfuncs;
	let record = ievaluate_args->record;
	let nrecord = ievaluate_args->nrecord;
	let popparam = &ievaluate_args->popparam;
	let logstream = ievaluate_args->logstream;

	char advan_memory[advanfuncs->advan_size];
	var advan = (ADVAN*)advan_memory;
	advanfuncs->construct(advan, advanfuncs);

	ADVAN_MODEL_MEMORY advanmem = { 0 };
	var predictvars = (PREDICTVARS*)advanmem._predictvars;
	var imodel = (IMODEL*)advanmem._imodel;
	assert((int)sizeof(advanmem._predictvars) >= advanfuncs->advanconfig->predictfields.size);
	assert((int)sizeof(advanmem._imodel) >= advanfuncs->advanconfig->imodelfields.size);

	/* iterate over individuals data */
	double lasttime = -DBL_MAX;
	let advanconfig = advanfuncs->advanconfig;
	let predict = advanconfig->predict;
	let predictall = advanconfig->predictall;
	let recordinfo = &advanfuncs->recordinfo;
	const RECORD* ptr = record;
	let id = RECORDINFO_ID(recordinfo, ptr);
	forcount(i, nrecord) {
		let evid = RECORDINFO_EVID(recordinfo, ptr);
		let time = RECORDINFO_TIME(recordinfo, ptr);

		let dv = RECORDINFO_DV(recordinfo, ptr);
		if (evid == 0 && !isfinite(dv))
			fatal(logstream, "DV is not finite for observation (EVID 0) at TIME (%.16e) for ID %f record %i\n", time, id, i);

		/* time must increase except for reset events */
		if (time < lasttime && evid != 3 && evid != 4)
			fatal(logstream, "TIME (%.16e) not monotonic for ID %f record %i, previous (%.16e)\n", time, id, i, lasttime);
		lasttime = time;

		/* compartment must be within the number of states */
		let cmt = RECORDINFO_CMT(recordinfo, ptr);
		if (cmt < 0 || cmt >= advanfuncs->nstate)
			fatal(logstream, "CMT (%.16e) not within number of states (%i) for ID %f record %i\n", cmt, id, i, advanfuncs->nstate);

		/* state should not be accessed outside of its limits */
		for (int j=advanfuncs->nstate; j<OPENPMX_STATE_MAX; j++)
			advan->state[j] = NAN;

		/* advance to the record time */
		let predictstate = advan_advance(advan, imodel, ptr, popparam);

		/* state always should remain finite */
		if (assure_state_finite(advan->state, advanfuncs->nstate) == 0) {
			forcount(j, advanfuncs->nstate)
				warning(logstream, "compartment %i state %.16e finite %i\n", j, advan->state[j], isfinite(advan->state[j]));
			fatal(logstream, "non-finite state for ID %f record %i time %f\n", id, i, time);
		}

		var yhat = 0.;
		if (evid == 0 || predictall)
			yhat = evaluate_yhat(imodel, &predictstate, popparam, predict, advanmem.errarray, predictvars);

		/* observations */
		if (evid == 0) {

			/* predictions for observations should be finite */
			if (isfinite(yhat) != 1)
				fatal(logstream, "non-finite YHAT for ID %f record %i\n", id, i);

			/* prediction variance of observations should be finite and positive */
			let yhatvar = evaluate_yhatvar(imodel, &predictstate, popparam, predict, advanmem.errarray, predictvars);
			if (isfinite(yhatvar) != 1)
				fatal(logstream, "non-finite YHATVAR for ID %f record %i\n", id, i);

			/* predictions with zero error are an error */
			if (evid == 0 && yhatvar == 0.)
				fatal(logstream, "zero YHATVAR for ID %f record %i\n", id, i);

			/* observation should be finite */
			if (isfinite(dv) != 1)
				fatal(logstream, "non-finite DV for ID %f record %i\n", id, i);

		/* dose or reset-and-dose event */
		} else if (evid == 1 || evid == 4) {

			let amt = RECORDINFO_AMT(recordinfo, ptr);
			if (!isfinite(amt))
				fatal(logstream, "ID %f record %i: AMT is not finite (%f)\n", id, i, amt);
			if (amt < 0.)
				fatal(logstream, "ID %f record %i: AMT less than zero\n", id, i);
			let rate = RECORDINFO_RATE(recordinfo, ptr);
			if (rate < 0.)
				fatal(logstream, "ID %f record %i: RATE less than zero\n", id, i);

			if (isfinite(amt) != 1 || amt == 0.)
				warning(logstream, "AMT is missing (%.16e) assumed 0 for dose event ID %f time %f record %i\n", amt, id, RECORDINFO_TIME(recordinfo, ptr), i);
			if (isfinite(rate) != 1)
				warning(logstream, "RATE is missing (%.16e) assumed 0 for dose event ID %f time %f record %i\n", rate, id, RECORDINFO_TIME(recordinfo, ptr), i);

			if (isfinite(dv) == 1 && dv != 0.)
				warning(logstream, "non-zero DV (%.16e) for dose event ID %f record %i\n", dv, id, i);

		/* reset event */
		} else if (evid == 3) {

			/* broken when AMT is NAN? */
			let amt = RECORDINFO_AMT(recordinfo, ptr);
			if (isfinite(amt) == 1 && amt != 0.)
				warning(logstream, "AMT is non-zero (%.16e) for reset event ID %f record %i\n", amt, id, i);

			let rate = RECORDINFO_RATE(recordinfo, ptr);
			if (isfinite(rate) == 1 && rate != 0.)
				warning(logstream, "RATE is non-zero (%.16e) for reset event ID %f record %i\n", rate, id, i);

			if (isfinite(dv) == 1 && dv != 0.)
				warning(logstream, "non-zero DV (%.16e) for reset event ID %f record %i\n", dv, id, i);

		/* other event type */
		} else {
			if (isfinite(dv) == 1 && dv != 0.)
				warning(logstream, "non-zero DV (%.16e) for non-observation for ID %f record %i\n", dv, id, i);
		}

		ptr = RECORDINFO_INDEX(recordinfo, ptr, 1);
	}
	advanfuncs->destruct(advan);
}

void individual_simulate(const IEVALUATE_ARGS* const ievaluate_args,
						 IMODEL* const imodel_saved,
						 PREDICTVARS* const predictvars_saved,
						 double* const istate,
						 double* const individ_yhat,
						 double* const individ_yhatvar,
						 double* const individ_pred,
						 const double* const isimerr)
{
	let advanfuncs = ievaluate_args->advanfuncs;
	let record = ievaluate_args->record;
	let nrecord = ievaluate_args->nrecord;

	char advan_memory[advanfuncs->advan_size];
	var advan = (ADVAN*)advan_memory;
	advanfuncs->construct(advan, advanfuncs);

	ADVAN_MODEL_MEMORY advanmem = { 0 }; /* errarray isnt used here, we use the one from the simulation */
	var predictvars = (PREDICTVARS*)advanmem._predictvars;
	var imodel = (IMODEL*)advanmem._imodel;

	/* iterate over individuals data */
	let advanconfig = advanfuncs->advanconfig;
	let predict = advanconfig->predict;
	let popparam = &ievaluate_args->popparam;
	let nsigma = popparam->nsigma;
	let recordinfo = &advanfuncs->recordinfo;
	const RECORD* ptr = record;
	forcount(i, nrecord) {
		let advanstate = advan_advance(advan, imodel, ptr, popparam);

		/* do prediction *without* the random error, put into yhat */
		/* yhat contains the prediction (*without* noise) */
		/* we dont to all predictions, otherwise it will be non-zero for infusions */
		var yhat = 0.;
		if (RECORDINFO_EVID(recordinfo, ptr) == 0)
			yhat = evaluate_yhat(imodel, &advanstate, popparam, predict, advanmem.errarray, predictvars);
		individ_yhat[i] = yhat;
		individ_yhatvar[i] = 0.;

		/* do prediction with the random error, put into pred */
		/* pred contains the prediction (*with* noise) */
		/* point to residual error with non-zero err values */
		let errarray = &isimerr[i * nsigma];
		var pred = 0.;
		if (RECORDINFO_EVID(recordinfo, ptr) == 0)
			pred = evaluate_yhat(imodel, &advanstate, popparam, predict, errarray, predictvars);
		individ_pred[i] = pred;

		/* save imodel, predictvars and state */
		let imodel_size = advanconfig->imodelfields.size;
		let imodelptr = (IMODEL*)((char*)imodel_saved + i * imodel_size);
		memcpy(imodelptr, imodel, imodel_size);

		let predictvars_size = advanconfig->predictfields.size;
		let predictvarsptr = (PREDICTVARS*)((char*)predictvars_saved + i * predictvars_size);
		memcpy(predictvarsptr, predictvars, predictvars_size);

		let nstate = advanfuncs->nstate;
		memcpy(&istate[i * nstate], advanstate.state, nstate * sizeof(double));

		ptr = RECORDINFO_INDEX(recordinfo, ptr, 1);
	}
	advanfuncs->destruct(advan);

}
