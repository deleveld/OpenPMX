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
#include "utils/c22.h"

/* NOTE: these functions must be thread safe on the level of an individual */

__attribute__ ((hot))
static inline double evaluate_yhat(const IMODEL* const imodel,
								   const PREDICTSTATE* const predictstate,
								   const POPPARAM* const popparam,
								   const IMODEL_PREDICT predict,
								   const double errarray[static OPENPMX_SIGMA_MAX],
								   PREDICTVARS* predictvars)
{
//	memset(predictvars, 0, OPENPMX_PREDICTVARS_MAX * sizeof(double));
		
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

			/* preserve zero errarray across function calls so we dont have
			 * to zero the entire errarray each time we call this */
			errarray[j] = 0;
			
			var deriv = (ya1 - ya2) / (above - below);
			yhatvar += (deriv*deriv) * sigma[j];
		}
	}
	return yhatvar;
}

static bool check_state(const double* const a, const int n, FILE* logstream, const bool _offset1)
{
	let record_offset = _offset1 ? 1 : 0;

	/* all state should be finite */
	int i;
	for (i=0; i<n; i++) {
		if (!isfinite(a[i])) {
			warning(logstream, "compartment %i state %f is not finite\n", i + record_offset, a[i]);
			return true;
		}
	}
	/* outside the used state should be NAN because we set that */
	for (i=n; i<OPENPMX_STATE_MAX; i++) {
		if (isfinite(a[i])) {
			warning(logstream, "compartment %i state %f should be NAN\n", i + record_offset, a[i]);
			return true;
		}
	}
	return false;
}

/* For how this is used see PAGE poster:
 * A comparison of methods for handling of data below the limit of quantification in NONMEM VI */
#include <gsl/gsl_cdf.h>
static inline double phi(const double x)
{
	return gsl_cdf_ugaussian_P(x);
}

typedef struct {
	double errarray[OPENPMX_SIGMA_MAX];
	double _predictvars[OPENPMX_PREDICTVARS_MAX];	/* will be cast to PREDICTVARS */
	double _imodel[OPENPMX_IMODEL_MAX];				/* will be cast to IMODEL */
} ADVAN_MODEL_MEMORY;

/* functions for Bae and Yim objective function Term 1 and Term 2 */
/* This is the core function evaluating an individual predictions to calculate
 * the objective function. Speeding up this function is quite important. */
 __attribute__ ((hot))
double individual_fasteval(const IEVALUATE_ARGS* const ievaluate_args)
{
	let advanfuncs = ievaluate_args->advanfuncs;
	let record = ievaluate_args->record;
	let nrecord = ievaluate_args->nrecord;
	let popparam = &ievaluate_args->popparam;

	char advan_memory[advanfuncs->advan_size];
	var advan = (ADVAN*)advan_memory;
	advanfuncs->construct(advan, advanfuncs);

	ADVAN_MODEL_MEMORY advanmem = { };
	var predictvars = (PREDICTVARS*)advanmem._predictvars;
	var imodel = (IMODEL*)advanmem._imodel;

	/* iterate over individuals data */
	let advanconfig = advanfuncs->advanconfig;
	let predict = advanconfig->predict;
	let recordinfo = &advanfuncs->recordinfo;
	let recordsize = recordinfo->dataconfig->recordfields.size;
	let no_dvlow_present = recordinfo->offsetDVLOW == -1;
	var obs_min2ll = 0.; /* separate sums to avoid loss of precision if magnitudes differ strongly */
	var obs_lndet = 0.;
	const RECORD* ptr = record;
	forcount(i, nrecord) {
		let predictstate = advan_advance(advan, imodel, ptr, popparam);
		if (RECORDINFO_EVID(recordinfo, ptr) == 0) {
			let yhat = evaluate_yhat(imodel, &predictstate, popparam, predict, advanmem.errarray, predictvars);
			let yhatvar = evaluate_yhatvar(imodel, &predictstate, popparam, predict, advanmem.errarray, predictvars);

			let dvlow = no_dvlow_present ? 0. : RECORDINFO_DVLOW(recordinfo, ptr);
			if (dvlow == 0.) {
				let dv = RECORDINFO_DV(recordinfo, ptr);
				let err = dv - yhat;
				obs_min2ll += (err * err) / yhatvar;
				obs_lndet += log(yhatvar);
			} else {
				let err = dvlow - yhat;
				obs_min2ll += -2. * log(phi(err / sqrt(yhatvar)));
				/* obs_lndet += 0.; This term does not take part I think. FIXME */
			}
		}
		ptr = RECORD_INDEX(ptr, recordsize, 1);	
	}
	advanfuncs->destruct(advan);

	return obs_min2ll + obs_lndet;
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
						 double* const ret_obs_lndet,
						 double* const ret_obs_min2ll)
{
	let advanfuncs = ievaluate_args->advanfuncs;
	let record = ievaluate_args->record;
	let nrecord = ievaluate_args->nrecord;
	let popparam = &ievaluate_args->popparam;

	char advan_memory[advanfuncs->advan_size];
	var advan = (ADVAN*)advan_memory;
	advanfuncs->construct(advan, advanfuncs);

	ADVAN_MODEL_MEMORY advanmem = { };
	var predictvars = (PREDICTVARS*)advanmem._predictvars;
	var imodel = (IMODEL*)advanmem._imodel;

	/* iterate over individuals data */
	let advanconfig = advanfuncs->advanconfig;
	let predict = advanconfig->predict;
	let predictall = advanconfig->predictall;
	let nstate = advanfuncs->nstate;
	let recordinfo = &advanfuncs->recordinfo;
	let no_dvlow_present = recordinfo->offsetDVLOW == -1;
	let imodel_size = advanconfig->imodelfields.size;
	let predictvars_size = advanconfig->predictfields.size;
	var obs_min2ll = 0.;
	var obs_lndet = 0.;
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

			let dv = RECORDINFO_DV(recordinfo, ptr);
			let dvlow = no_dvlow_present ? 0. : RECORDINFO_DVLOW(recordinfo, ptr);
			if (dvlow == 0.) {
				let err = dv - yhat;
				obs_min2ll += (err * err) / yhatvar;
				obs_lndet += log(yhatvar);
			} else {
				let err = dvlow - yhat;
				obs_min2ll += -2. * log(phi(err / sqrt(yhatvar)));
			}
		} else {
			/* yhatvar must be set to zero for non-observations since these
			 * will be used to calculate the individual covariance matrix */
			if (YHATVAR)
				YHATVAR[i] = 0.;
		}
		ptr = RECORDINFO_INDEX(recordinfo, ptr, 1);
	}
	advanfuncs->destruct(advan);

	if (ret_obs_min2ll)
		*ret_obs_min2ll = obs_min2ll;
	if (ret_obs_lndet) 
		*ret_obs_lndet = obs_lndet;
}

void individual_checkout(const IEVALUATE_ARGS* const ievaluate_args)
{
	let advanfuncs = ievaluate_args->advanfuncs;
	let record = ievaluate_args->record;
	let nrecord = ievaluate_args->nrecord;
	let popparam = &ievaluate_args->popparam;
	let logstream = ievaluate_args->logstream;

	/* if RATE exists but AMT does not it is an error */
	if (advanfuncs->recordinfo.offsetRATE != -1 && advanfuncs->recordinfo.offsetAMT == -1) 
		fatal(0, "error: RATE exists but AMT does not\n");

	char advan_memory[advanfuncs->advan_size + 1000];
	memset(advan_memory, 0, advanfuncs->advan_size + 1000);
	var advan = (ADVAN*)advan_memory;
	advanfuncs->construct(advan, advanfuncs);

	ADVAN_MODEL_MEMORY advanmem = { };
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
	let no_dvlow_present = recordinfo->offsetDVLOW == -1;
	const RECORD* ptr = record;
	let id = RECORDINFO_ID(recordinfo, ptr);
	if (id != floor(id))
		warning(logstream, "ID (%f) should probably be an integer\n", id);
	let _offset1 = recordinfo->_offset1;
	let record_offset = _offset1 ? 1 : 0;
	forcount(i, nrecord) {
		let time = RECORDINFO_TIME(recordinfo, ptr);
		let dv = RECORDINFO_DV(recordinfo, ptr);
		let evid = RECORDINFO_EVID(recordinfo, ptr);

		let cmt = RECORDINFO_CMT(recordinfo, ptr);
		if (cmt != floor(cmt))
			warning(logstream, "CMT (%f) should probably be an integer: ID %f time %f record %i\n", cmt, id, time, i + record_offset);

		/* time must increase except for reset events */
		if (time < lasttime && evid != 3 && evid != 4)
			fatal(logstream, "TIME not monotonic: ID %f time %f record %i previous %f\n", id, time, i + record_offset, lasttime);
		lasttime = time;

		/* state should not be accessed outside of its limits */
		for (int j=advanfuncs->nstate; j<OPENPMX_STATE_MAX; j++)
			advan->state[j] = NAN;

		/* check record before advance */
		/* observations */
		if (evid == 0) {
			if (dv == 0.)
				warning(logstream, "DV zero for observation: ID %f time %f record %i\n", id, time, i + record_offset);

			/* observation should be finite */
			if (!isfinite(dv))
				fatal(logstream, "DV not finite for observation: ID %f time %f record %i\n", id, time, i + record_offset);

			let amt = RECORDINFO_AMT(recordinfo, ptr);
			if (isfinite(amt) && amt != 0.)
				warning(logstream, "AMT non-zero (%f) for observation: ID %f time %f record %i\n", amt, id, time, i + record_offset);

			let rate = RECORDINFO_RATE(recordinfo, ptr);
			if (isfinite(rate) && rate != 0.)
				warning(logstream, "RATE non-zero (%f) for observation: ID %f time %f record %i\n", rate, id, time, i + record_offset);

		/* dose or reset-and-dose event */
		} else if (evid == 1 || evid == 4) {

			/* compartment must be within the number of states */
			let cmt = RECORDINFO_CMT_0offset(recordinfo, ptr);
			if (cmt < 0 || cmt >= advanfuncs->nstate) {
				let _cmt = RECORDINFO_CMT(recordinfo, ptr);
				fatal(logstream, "CMT (%i) not within number of states (%i): ID %f time %f record %i\n", _cmt, advanfuncs->nstate, id, time, i + record_offset);
			}

			let amt = RECORDINFO_AMT(recordinfo, ptr);
			if (!isfinite(amt))
				fatal(logstream, "AMT not finite (%f): ID %f time %f record %i\n", amt, id, time, i + record_offset);
			if (amt <= 0.)
				fatal(logstream, "AMT less than or equal to zero (%f): ID %f time %f record %i \n", amt, id, time, i + record_offset);
			let rate = RECORDINFO_RATE(recordinfo, ptr);
			if (rate < 0.)
				fatal(logstream, "RATE less than zero (%f): ID %f record %i\n", rate, id, i + record_offset);

			if (!isfinite(amt)|| amt == 0.)
				warning(logstream, "AMT missing (%f) assumed 0 for dose: ID %f time %f record %i\n", amt, id, time, i + record_offset);
			if (!isfinite(rate))
				warning(logstream, "RATE missing (%f) assumed 0 for dose: ID %f time %f record %i\n", rate, id, time, i + record_offset);

			if (isfinite(dv) == 1 && dv != 0.)
				warning(logstream, "DV non-zero (%f) for dose: ID %f time %f record %i\n", dv, id, time, i + record_offset);

		/* reset event */
		} else if (evid == 3) {

			let amt = RECORDINFO_AMT(recordinfo, ptr);
			if (isfinite(amt) && amt != 0.)
				warning(logstream, "AMT non-zero (%f) for reset: ID %f time %f record %i\n", amt, id, time, i + record_offset);

			let rate = RECORDINFO_RATE(recordinfo, ptr);
			if (isfinite(rate) && rate != 0.)
				warning(logstream, "RATE non-zero (%f) for reset: ID %f time %f record %i\n", rate, id, time, i + record_offset);

			if (isfinite(dv) && dv != 0.)
				warning(logstream, "DV non-zero (%f) for reset: ID %f time %f record %i\n", dv, id, time, i + record_offset);

		/* other event type */
		} else {
			if (!isfinite(dv) && dv != 0.)
				warning(logstream, "non-zero DV (%f) for non-observation: ID %f time %f record %i\n", dv, id, time, i + record_offset);
		}

		/* dvlow implies observation */
		if (evid != 0) {
			let dvlow = no_dvlow_present ? 0. : RECORDINFO_DVLOW(recordinfo, ptr);
			if (dvlow != 0.) 
				fatal(logstream, "DVLOW (%f) present for non-observation: ID %f record %i\n", dvlow, id, i + record_offset);
		}

		/* state always should remain finite */
		if (check_state(advan->state, advanfuncs->nstate, logstream, _offset1)) 
			fatal(logstream, "non-finite state before advance: ID %f time %f record %i\n", id, time, i + record_offset);

		/* advance to the record time */
		let predictstate = advan_advance(advan, imodel, ptr, popparam);

		/* now check record after advance */

		/* state always should remain finite */
		if (check_state(advan->state, advanfuncs->nstate, logstream, _offset1)) 
			fatal(logstream, "non-finite state after advance: ID %f time %f record %i\n", id, time, i + record_offset);

		/* predictions should be finite */
		var yhat = 0.;
		if (evid == 0 || predictall) {
			yhat = evaluate_yhat(imodel, &predictstate, popparam, predict, advanmem.errarray, predictvars);
			if (!isfinite(yhat))
				fatal(logstream, "YHAT non-finite: ID %f time %f record %i\n", id, time, i + record_offset);
		}

		/* observations */
		if (evid == 0) {

			/* prediction variance of observations should be finite and positive */
			let yhatvar = evaluate_yhatvar(imodel, &predictstate, popparam, predict, advanmem.errarray, predictvars);
			if (!isfinite(yhatvar))
				fatal(logstream, "YHATVAR non-finite: ID %f time %f record %i\n", id, time, i + record_offset);

			/* predictions with zero error are an error */
			if (evid == 0 && yhatvar == 0.)
				fatal(logstream, "YHATVAR zero (%f): ID %f time %f record %i\n", yhatvar, id, time, i + record_offset);
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

	ADVAN_MODEL_MEMORY advanmem = { }; /* errarray isnt used here, we use the one from the simulation */
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

IEVALUATE_ARGS ievaluate_args_init(const RECORD* const record,
								   const int nrecord,
								   const ADVANFUNCS* const advanfuncs,
								   const double* const theta,
								   const int ntheta,
								   const double* const eta,
								   const int nomega,
								   const double* const sigma,
								   const int nsigma,
								   FILE* logstream)
{	
	return (IEVALUATE_ARGS) {
		.record = record,
		.nrecord = nrecord,
		.advanfuncs = advanfuncs,
		.popparam = { 
			.theta = theta,
			.ntheta = ntheta,
			.eta = eta,
			.nomega = nomega,
			.sigma = sigma,
			.nsigma = nsigma,
			.nstate = advanfuncs->nstate 
		},
		.logstream = logstream,
	};
}
