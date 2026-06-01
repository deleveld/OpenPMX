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

/// This file implements at Target-Controlled-Infusion dose controller.

#include <string.h>
#include <values.h>
#include <assert.h>
 
#include "advan/advan.h"
#include "utils/c22.h"

/// Jona Joachim (jona@joachim.cc) Contributed the code in `src/advan/ministan` which is 
/// a modernized version of [Stanpump](https://opentci.org/code/stanpump).
/* include the ministan source directly */ 
#include "advan/ministan/common.h"
#include "advan/ministan/cube.c"
#include "advan/ministan/udfs.c"
#include "advan/ministan/helpers.c"
#include "advan/ministan/virtual_model.c"
#include "advan/ministan/model.c"
#include "advan/ministan/find_peak.c"
#include "advan/ministan/cfg.c"

typedef struct TCICONTROL {
	Config cfg;
	const double max_rate;
	const double peak_time;
	const int cmt_0;	/* 0-based compartment */
	double next_time;	/* in TIME units */
	double totamt; 		/* cumulative dose */
	ADVANINFUSION last_infusion;
	const TCICONFIG last_tciconfig;
} TCICONTROL;

bool pmx_advan_tci_started(const ADVANSTATE* advanstate)
{
	var tcicontrol = advanstate->advan->tcicontrol;
	return tcicontrol ? true : false;
}

void pmx_advan_tci_init(const ADVANSTATE* advanstate, const TCICONFIG* const tciconfig)
{
	/* do nothing if already initialzed, and the same settings */
	var tcicontrol = advanstate->advan->tcicontrol;
	if (tcicontrol) {

		/* its probably an error to change the config during TCI */
		if (memcmp(&tcicontrol->last_tciconfig,
				   tciconfig,
				   sizeof(TCICONFIG)) != 0) {
			printf("fatal: %s: TCI config does not match current\n", __func__); 
			exit(EXIT_FAILURE);
		}
		return;
	}
		
	/* if cmt=0 for 1-offset then we correct to first (0) compartment */
	var cmt_0 = (advanstate->advan->advanfuncs->recordinfo.dataconfig->_offset1) ? tciconfig->cmt - 1 : tciconfig->cmt;
	if (cmt_0 < 0)
		cmt_0 = 0;

	/* Config is binary copyable so we can construct it in place and 
	 * copy it to malloc memory like this */
	var tcicfg = (TCICONTROL) {
		.cfg = {
			.k10 = tciconfig->k10,
			.k12 = tciconfig->k12,
			.k13 = tciconfig->k13,
			.k21 = tciconfig->k21,
			.k31 = tciconfig->k31,
			.ke0 = tciconfig->ke0,
			.vc  = tciconfig->vc,
			.target_effect = tciconfig->target_effect,
			.delta_seconds = 10,
		},
		/* .peak_time set after cfg made */
		.cmt_0 = cmt_0,
		.max_rate = tciconfig->max_rate,
		.next_time = DBL_MAX,
		.totamt = 0.,
		.last_tciconfig = *tciconfig,
	};
	cfg_init(&tcicfg.cfg);
	
	tcicontrol = mallocvar(TCICONTROL, 1);
	memcpy(tcicontrol, &tcicfg, sizeof(TCICONTROL));
	
	/* write back into advan, so we access it in other functions */
	advanstate->advan->tcicontrol = tcicontrol;
}

static inline void add_infusion(ADVANINFUSION *arr, int *size, const ADVANINFUSION* const newinf)
{
	arr[*size] = *newinf;
    (*size)++;
    assert(*size < OPENPMX_SIMULINFUSION_MAX);
}

double pmx_advan_tci_target(const ADVANSTATE* advanstate, const double target)
{
	var tcicontrol = advanstate->advan->tcicontrol;
	if (!tcicontrol) {
		fprintf(stderr, "fatal: %s: tci not initialized\n", __func__);
		exit(EXIT_FAILURE);
	}

	/* setting a target < 0 turn off the TCI controller */
	let totamt = tcicontrol->totamt;
	if (target < 0) {
		free(tcicontrol);
		advanstate->advan->tcicontrol = 0;
		return totamt;
	}

	let delta_seconds = tcicontrol->cfg.delta_seconds;
	let duration = delta_seconds / 60.;

	/* we only have to do something when the previous infusion has stopped */
	let now = advanstate->advan->time;
	var next_time = tcicontrol->next_time;
	let last = &tcicontrol->last_infusion;
	if (next_time > now && next_time != DBL_MAX) {
		let elapsed = fmin(now - last->start, duration);
		return tcicontrol->totamt + last->rate * elapsed;
	}

	/* the previous infusion has completed, accumulate its full amount */
	if (next_time != DBL_MAX)
		tcicontrol->totamt += last->amt;
	tcicontrol->last_infusion = (ADVANINFUSION){ };
		
	/* Use stanpump to calculate the rate */
	var cfg = &tcicontrol->cfg;
	var rate_s = calculate_rate(cfg, target);

	/* apply maximum rate */
	if (tcicontrol->max_rate > 0.)
		rate_s = fmin(rate_s, tcicontrol->max_rate);

	/* use the right rate */
	advance_rate(cfg, rate_s);

	/* apply the rate as an infusion, starting now */
	/* make infusion a tiny bit faster to ensure correct ordering
	 * TODO: This should probably be made more elegant. Im not sure if
	 * we would otherwise have to then sort the infusions, or
	 * prefer to handle certian infusions before records etc */
	next_time = now + duration - 1e-12;
	if (rate_s != 0.) {
		let amt = rate_s * delta_seconds;
		let rate = amt / duration;
		let advan = advanstate->advan;
		let inf = (ADVANINFUSION) {
			.cmt = tcicontrol->cmt_0,
			.amt = amt,
			.rate = rate,
			.start = now,
			.end = next_time,
		};
		add_infusion(advan->infusions, &advan->ninfusions, &inf);
		tcicontrol->last_infusion = inf;
	}

	/* come back later to continue the TCI and recalculate the rate */
	tcicontrol->next_time = next_time;
	pmx_advan_inittime(advanstate, next_time);

	return tcicontrol->totamt;
}

double pmx_advan_tci_plasma_conc(const ADVANSTATE* advanstate)
{
	var tcicontrol = advanstate->advan->tcicontrol;
	if (!tcicontrol) {
		fprintf(stderr, "fatal: %s: tci not initialized\n", __func__);
		exit(EXIT_FAILURE);
	}
	return tcicontrol->cfg.plasma_conc;
}

double pmx_advan_tci_effect_conc(const ADVANSTATE* advanstate)
{
	var tcicontrol = advanstate->advan->tcicontrol;
	if (!tcicontrol) {
		fprintf(stderr, "fatal: %s: tci not initialized\n", __func__);
		exit(EXIT_FAILURE);
	}
	return tcicontrol->cfg.effect_conc;
}


