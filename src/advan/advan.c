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

#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "advan.h"
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
//	vector_reserve(advan->infusions, OPENPMX_SIMULINFUSION_MAX);
	vector_init_buffer(advan->infusions, advan->_infusions_buffer, OPENPMX_SIMULINFUSION_MAX);
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

	/* reset the state on first start or on EVID==3 or EVID==4 events */
	let evid = RECORDINFO_EVID(recordinfo, record);
	let reset_state = (advan->initcount == 0 || evid == 3 || evid == 4);

	if (reset_state || advanconfig->firstonly == false) {
		if (reset_state) {
			/*	If we dont zero the model it keeps the values from the previous
			 * time we called. This seems to be somewhat elegant in the user code.
			 * Since the model retains values from the previous record. 
			 * memset(imodel, 0, advanfuncs->advanconfig->imodelfields.size); */
 
			memset(state, 0, OPENPMX_STATE_MAX * sizeof(double));
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
							.statetime = advan->time,
							.record = record,
							.state = state },
						  popparam);
		++advan->initcount;
	}

	/* now handle any doses, if this is a dose add it to the list */
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
			let v = &advan->infusions.data[i];
			if (v->end <= currenttime && v->rate > 0.) {
				vector_remove(advan->infusions, i, 1);
				--i;
			}
		}

		/* find the first place we have to stop at going to the next record */
		let recordtime = RECORDINFO_TIME(recordinfo, record);
		double intervalstop = recordtime;
		forvector(i, advan->infusions) {
			let v = &advan->infusions.data[i];
			if (v->end < intervalstop)
				intervalstop = v->end;
			if (v->start > currenttime && v->start < intervalstop)
				intervalstop = v->start;
		}

		/* do we advance time at all */
		if (intervalstop > currenttime) {

			/* collect the infusion rates between now and the intervalstop */
			double totalrates[OPENPMX_STATE_MAX] = { 0 };
			forvector(i, advan->infusions) {
				let v = &advan->infusions.data[i];
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
		/* a fake boilus with cmt == -1 means an extra call to init */
		bool need_reset_now = false;
		bool need_init_now = false;
		forvector(i, advan->infusions) {
			let v = &advan->infusions.data[i];
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
								.statetime = advan->time,
								.record = record,
								.state = state,
							  },
							  popparam);
		}

	/* keep going till we have the state equal to the required record time */
	} while (advan->time < RECORDINFO_TIME(recordinfo, record));

	/* return the state useful for calling predict */
	return (PREDICTSTATE) {
		.state = state,
		.record = record,
	};
}

int pmx_advan_initcount(const ADVANSTATE* const advanstate)
{
	var advan = advanstate->advan;
	return advan->initcount;
}

void pmx_advan_amtlag(const ADVANSTATE* const advanstate, const int cmt, const double t)
{
	var advan = advanstate->advan;
	let advanfuncs = advan->advanfuncs;
	let nstate = advanfuncs->nstate;

	assert(cmt < nstate);
	advan->amtlag[cmt] = t;
}

void pmx_advan_bioaval(const ADVANSTATE* const advanstate, const int cmt, const double f)
{
	var advan = advanstate->advan;
	let advanfuncs = advan->advanfuncs;
	let nstate = advanfuncs->nstate;

	assert(cmt < nstate);
	advan->bioavail[cmt] = f;
}

void pmx_advan_inittime(const ADVANSTATE* const advanstate, const double t)
{
	var advan = advanstate->advan;

	/* ignores setting in the past */
	if (t <= advan->time)
		return;

	/* if its already on the list then ignore */
	bool found_already = false;
	forvector(i, advan->infusions) {
		let v = &advan->infusions.data[i];
		if (v->cmt == -1. && v->start == t && v->end == t)
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

void pmx_advan_state_init(const ADVANSTATE* const advanstate, const int cmt, const double v)
{
	var advan = advanstate->advan;
	if (advan->initcount == 0) {
		printf("advan initcount %i time %f\n", advan->initcount, advan->time);
		let advanfuncs = advan->advanfuncs;
		let nstate = advanfuncs->nstate;
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

