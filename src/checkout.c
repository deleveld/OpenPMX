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

#include "checkout.h"
#include "ievaluate.h"
#include "scatter.h"
#include "print.h"
#include "utils/c22.h"
#include "utils/various.h"

#include "openpmx_compile_options.h"

/* NOTE: this function must be thread safe on the level of an individual */
static void idata_checkout_thread(INDIVID* const individ,
								  const ADVANFUNCS* const advanfuncs,
								  const POPMODEL* const popmodel,
								  const NONZERO* const nonzero,
								  const OPTIONS* const options,
								  const SCATTEROPTIONS* const scatteroptions)
{
	(void)nonzero;
	(void)options;

	struct timespec t3;
	clock_gettime(CLOCK_REALTIME, &t3);

	/* checkout is at ETA at 0 */
	double eta[OPENPMX_OMEGA_MAX] = { 0 };
	let ievaluate_args = (IEVALUATE_ARGS) {
		.record = individ->record,
		.nrecord = individ->nrecord,
		.advanfuncs = advanfuncs,
		.popparam = popparam_init(popmodel, advanfuncs, eta),
		.logstream = scatteroptions->logstream,
	};
	individual_checkout(&ievaluate_args);

	individ->ineval += 1;
	individ->neval += 1;
	timespec_duration(&t3, &individ->eval_msec);
}

void idata_checkout(const IDATA* const idata,
					const ADVANFUNCS* const advanfuncs,
					const POPMODEL* const popmodel,
					const OPTIONS* const options,
					FILE* logstream)
{
	info(logstream, "checkout begin\n");

	SCATTEROPTIONS scatteroptions = { 0 };
	scatteroptions.checkout_errors = true;
	scatteroptions.logstream = logstream;
	scatter_threads(idata, advanfuncs, popmodel, 0, options, &scatteroptions, idata_checkout_thread);

	info(logstream, "checkout end\n");

/* checkout is only called just before a run to we dont have to stop the
 * threads as they will be used shortly anyway. */
/*#ifdef OPENPMX_PARALLEL_PTHREADS
	pthreads_cleanup(logstream);
#endif */
}

