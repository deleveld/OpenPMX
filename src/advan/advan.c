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
 
/// The user creates an ADVANFUNCS object via a method named something 
/// like: 
/// `ADVANFUNCS* pmx_advan_FOO(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig)`
/// This function allocates and constructs an object containing function 
/// pointers to construct, advance across records, and destruct an 
/// advancer object. This is used in ievaluate.c.
	
/// This file is the base class ADVAN which contains information that 
/// all of the advancer objects have in common, for example: time, 
/// state, lag, bioavailability, and the infusions running at that 
/// moment. Probably the most important function is `advan_advance()` 
/// which handles the infusion AMT, RATE, EVID and the various special 
/// flags like CMT. 

#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "print.h"
#include "advan/advan.h"
#include "utils/c22.h"

void advan_base_construct(ADVAN* advan, const ADVANFUNCS* advanfuncs)
{
	memset(advan, 0, advanfuncs->advan_size); /* zero the full size */

	advan->advanfuncs = advanfuncs;
	advan->time = NAN;
	advan->initcount = 0;

	for (int i=0; i<OPENPMX_STATE_MAX; i++)
		advan->bioavail[i] = 1.;

	/* this allows running without any allocations. TODO: There should be a warning
	 * if more there an unexpected number of overlapping infusions. Right now vector_reserve()
	 * would assert() */
	vector_reserve(advan->infusions, OPENPMX_SIMULINFUSION_MAX);
}

void advan_base_destruct(ADVAN* advan)
{
	vector_free(advan->infusions);
}

ADVANFUNCS* advanfuncs_alloc(const DATACONFIG* dataconfig, const ADVANCONFIG* const advanconfig)
{
	assert(advanconfig->method);
	return advanconfig->method(dataconfig, advanconfig); /* This has to construct the RECORDINFO */
}

void advanfuncs_free(ADVANFUNCS* advanfuncs)
{
	free(advanfuncs);
}

PREDICTSTATE advan_advance(ADVAN* const advan,
						   IMODEL* const imodel,
						   const RECORD* const record,
						   const POPPARAM* const popparam)
{
	let advanfuncs = advan->advanfuncs;
	let advanconfig = advanfuncs->advanconfig;
	let recordinfo = &advanfuncs->recordinfo;
	var state = advan->state;
	let nstate = advanfuncs->nstate;

/// ### Advancing through an individuals records

/// While advancing the init() function is called (which calls the users
/// IMODEL() code in openpmxtran) when:
///
/// + First record of an individual
/// + EVID is 3 (reset event) or EVID is 4 (reset-and-dose event)
	/* reset the state on first start or on EVID==3 or EVID==4 events */
	let evid = RECORDINFO_EVID(recordinfo, record);
	let reset_state = (advan->initcount == 0 || evid == 3 || evid == 4);

	if (reset_state || advanconfig->firstonly == false) {
		if (reset_state) {
			/*	If we dont zero the model it keeps the values from the previous
			 * time we called. This seems to be somewhat elegant in the user code.
			 * Since the model retains values from the previous record. 
			 * memset(imodel, 0, advanfuncs->advanconfig->imodelfields.size); */
			if (nstate)
				memset(state, 0, nstate * sizeof(double));
			advan->time = RECORDINFO_TIME(recordinfo, record);
			vector_resize(advan->infusions, 0);

			/* for reset and reset-and-dose we should treat it like a discontinuitiy
			 * This helps the ODE solver a lot I expect */
			if (advan->initcount != 0 && advanfuncs->reset)
				advanfuncs->reset(advan, 1);
		}
		advanconfig->init(imodel,
						  &(ADVANSTATE) {
							.advan = advan,
							.current = {
								.statetime = advan->time,
								.record = record,
								.state = state,
								.popparam = popparam,
							},
						  });
		++advan->initcount;
	}

/// + For dose records EVID is 1 or EVID is 4 
/// + At the time given by a call to `pmx_advan_inittime()` which is 
/// `INITTIME()` in openpmxtran
///
/// The moment a dose record is encountered it is saved in a buffer with
/// the time of the RECORD and the lag valid at the moment of 
/// processing. Subsequent changes in lag do not change the moment when 
/// the dose is actually applied.
	if (evid == 1 || evid == 4) {
		let amt = RECORDINFO_AMT(recordinfo, record);
		if (amt > 0.) {
			let cmt = RECORDINFO_CMT_0offset(recordinfo, record);
			let lagtime = advan->amtlag[cmt];
			let start = RECORDINFO_TIME(recordinfo, record) + lagtime;
			let rate = RECORDINFO_RATE(recordinfo, record);
			let duration = (rate != 0.) ? (amt / rate) : (0.);
			let end = start + duration;

			vector_append(advan->infusions,
						  (ADVANINFUSION) {
							.cmt = cmt,
							.amt = amt,
							.rate = rate,
							.start = start,
							.end = end } );
		}
	}
	
	/* advance through time and handle infusions and doses when they start
	 * or stop up until the time of the current record */
	assert(advan->time <= RECORDINFO_TIME(recordinfo, record));
	do {
		let currenttime = advan->time;

		/* remove infusions that are in the past */
		forvector(i, advan->infusions) {
			let v = &advan->infusions.ptr[i];
			if (v->end <= currenttime && v->rate > 0.) {
				vector_remove(advan->infusions, i, 1);
				--i;
			}
		}

/// Advancing is done in "steps" to the earliest of:
///
/// + The next RECORD time is achieved
/// + A previous infusion starts or stopped
/// + A bolus dose is given
/// + If `pmx_advan_inittime()` has been called, which is INITTIME() in 
/// openpmxtran
		/* find the first place we have to stop at going to the next record */
		let recordtime = RECORDINFO_TIME(recordinfo, record);
		double intervalstop = recordtime;
		forvector(i, advan->infusions) {
			let v = &advan->infusions.ptr[i];
			if (v->end < intervalstop)
				intervalstop = v->end;
			if (v->start > currenttime && v->start < intervalstop)
				intervalstop = v->start;
		}

		/* do we advance time at all */
		if (intervalstop > currenttime) {

			/* collect the infusion rates between now and the intervalstop */
			double totalrates[OPENPMX_STATE_MAX] = { };
			forvector(i, advan->infusions) {
				let v = &advan->infusions.ptr[i];
				if (v->start <= currenttime && v->end >= currenttime && v->rate > 0.) {
					let cmt = v->cmt;
					let rate = v->rate;
					var bioavail = advan->bioavail[cmt];
					totalrates[cmt] += rate * bioavail;
				}
			}
			advanfuncs->interval(advan, imodel, record, state, popparam, intervalstop, totalrates);
			advan->time = intervalstop;	/* force and update to time, although the advancer_object might do this as well */
		}

		/* are any bolus needed to be given right now */
		/* we take them off the infusion list so they wont be seen again */
		/* a fake bolus with cmt == -1 means an extra call to init */
		bool need_reset_now = false;
		bool need_init_now = false;
		forvector(i, advan->infusions) {
			let v = &advan->infusions.ptr[i];
			if (v->start == currenttime && v->rate == 0.) {
				assert(v->start == v->end);
				let record_cmt = v->cmt;
				if (v->cmt == -1.) {		/* This means call init again */
					need_init_now = true;
				} else {
					let record_amt = v->amt;
					var bioavail = advan->bioavail[record_cmt];
					state[record_cmt] += record_amt * bioavail;
				}
				vector_remove(advan->infusions, i, 1);
				--i;
				need_reset_now = true;
			}
		}
		if (need_reset_now && advanfuncs->reset)
			advanfuncs->reset(advan, 0);

		/* extra call to init */
		if (need_init_now) {
			advanconfig->init(imodel,
							&(ADVANSTATE) {
								.advan = advan,
								.current = {
									.statetime = advan->time,
									.record = record,
									.state = state,
									.popparam = popparam,
								},
							});
		}

	/* keep going till we have the state equal to the required record time */
	} while (advan->time < RECORDINFO_TIME(recordinfo, record));

	/* return the state useful for calling predict */
	return (PREDICTSTATE) {
		.statetime = advan->time,
		.state = state,
		.record = record,
		.popparam = popparam,
	};
}

/* We dont have to worry that this function might be called outside of
 * init function because the ADVANSTATE is only available there */
void pmx_advan_amtlag(const ADVANSTATE* advanstate, const int cmt, const double t)
{
	var advan = advanstate->advan;
	let advanfuncs = advan->advanfuncs;
	let nstate = advanfuncs->nstate;
	
/// Calling `pmx_advan_amtlag()` (in openpmxtran this is `ALAG()`) 
/// delays the application of subsequent dosing. It does not change the 
/// moment of application of doses that are already in the buffer.
	assert(!isnan(t));
	assert(cmt >= 0);
	assert(cmt < nstate);
	
	advan->amtlag[cmt] = t;
}

/* We dont have to worry that this function might be called outside of
 * init function because the ADVANSTATE is only available there */
void pmx_advan_bioaval(const ADVANSTATE* advanstate, const int cmt, const double f)
{
	var advan = advanstate->advan;
	let advanfuncs = advan->advanfuncs;
	let nstate = advanfuncs->nstate;

	assert(!isnan(f));
	assert(cmt >= 0);
	assert(cmt < nstate);
	
	advan->bioavail[cmt] = f;
}

/* We dont have to worry that this function might be called outside of
 * init function because the ADVANSTATE is only available there */
void pmx_advan_inittime(const ADVANSTATE* advanstate, const double t)
{
	var advan = advanstate->advan;

	assert(!isnan(t));

	/* ignores setting in the past */
	if (t <= advan->time)
		return;

	/* its kind of not clear whether we should or shouldnt ignore the extra
	 * inittime if another kind of event occurs at the intended time, for
	 * example an infusion */
	bool found_already = false;
	forvector(i, advan->infusions) {
		let v = &advan->infusions.ptr[i];
		if (v->start == t)
			found_already = true;
	}

	if (!found_already) {
		vector_append(advan->infusions,
					  (ADVANINFUSION) {
						.cmt = -1,
						.amt = 0.,
						.rate = 0.,
						.start = t,
						.end = t });
	}
}

void pmx_advan_state_init(const ADVANSTATE* advanstate, const int cmt, const double v)
{
	var advan = advanstate->advan;
	let advanfuncs = advan->advanfuncs;
	let nstate = advanfuncs->nstate;

	assert(!isnan(v));
	assert(cmt >= 0);
	assert(cmt < nstate);

	if (advan->initcount == 0) {
		assert(cmt < nstate);
		advan->state[cmt] = v;
	}
}

/*
double record_variable(const char* name, const POPPARAM* const popparam, const RECORD* record, int* offset)
{
	if (offset && *offset >= 0)
		return DATA_FIELD(record, *offset);

	let realoffset = structinfo_find_offset(name, popparam->recordfields);
	if (offset)
		*offset = realoffset;

	assert(realoffset >= 0);

	return DATA_FIELD(record, realoffset);
}
*/

