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
#include <math.h>

#include "advan.h"
#include "utils/c22.h"

typedef struct {
	ADVANFUNCS advanfuncs;
} ADVANTABLE_WRAPPER;

static void advancer_wrapper_info(const struct ADVANFUNCS* const advanfuncs, FILE* f)
{
	let wrapper = (const ADVANTABLE_WRAPPER*)advanfuncs;
	let inner = advanfuncs->advanconfig->args.inner.method;
	
	fprintf(f, "advan wrapper\n");
	
	fprintf(f, "advan stepper: %s\n", libgsl->steptype_name);
	fprintf(f, "advan nstate: %i\n", advanfuncs->nstate);
	fprintf(f, "advan abstol: %g\n", libgsl->abstol);
	fprintf(f, "advan reltol: %g\n", libgsl->reltol);
	fprintf(f, "advan hstart: %g\n", libgsl->hstart);


	
}

static void advancer_wrapper_construct(ADVAN* advan, const ADVANFUNCS* const advanfuncs)
{
}

static void advancer_wrapper_destruct(ADVAN * advan)
{
}

static void advancer_wrapper_reset(ADVAN * advan, const int full)
{
	(void) advan;
	(void) full;
}

static void advancer_wrapper_advance_interval(ADVAN* advan,
											  const IMODEL* const imodel,
											  const RECORD* const record,
											  double* const state,
											  const POPPARAM* const popparam,
											  const double endtime,
											  const double* rates)
{
}

ADVANFUNCS* pmx_advan_wrapper(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig)
{
	let inner = advanconfig->args.inner.method;
	var ret = advanconfig->args.inner.method(dataconfig, advanconfig);

	/* make extra room, this will definitly be enough */
	let new_advan_size = ret->size + sizeof(ADVANTABLE_WRAPPER);

	ret = realloc(ret, new_advan_size);
	ret->size = ret->size

	ret->info = advancer_wrapper_info;
	ret->construct = advancer_wrapper_construct;
	ret->destruct = advancer_wrapper_destruct;
	ret->reset = advancer_wrapper_reset;
	ret->interval = advancer_wrapper_interval;
	
	return ret;
}

