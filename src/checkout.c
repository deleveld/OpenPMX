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

/// This file implements the high-level checkout functions. The actual
/// checkout function is individual_checkout() which is implemented in
/// ievaluate.c. It advances the individual records with ETA zet to 0
/// and has lots of checks to detect errors. This is important because
/// the other advancing of the individual records for objective function
/// calculation, prediction, etc do not have any check to maximize
/// speed.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "checkout.h"
#include "ievaluate.h"
#include "scatter.h"
#include "print.h"
#include "advan/advan.h"
#include "utils/c22.h"
#include "utils/various.h"

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
	clock_gettime(CLOCK_MONOTONIC, &t3);

	/* checkout is at ETA at 0 */
	double eta[OPENPMX_OMEGA_MAX];
	forcount(i, OPENPMX_OMEGA_MAX)
		eta[i] = NAN;
	forcount(i, popmodel->nomega)
		eta[i] = 0.;
	
	let ievaluate_args = ievaluate_args_init(individ->record,
											 individ->nrecord,
											 advanfuncs,
											 popmodel->theta,
											 popmodel->ntheta,
											 eta,
											 popmodel->nomega,
											 popmodel->sigma,
											 popmodel->nsigma,
											 scatteroptions->logstream);
	individual_checkout(&ievaluate_args);

	individ->ineval += 1;
	timespec_duration(&t3, &individ->eval_msec);
}

void idata_checkout(const IDATA* const idata,
					const ADVANFUNCS* const advanfuncs,
					const POPMODEL* const popmodel,
					const OPTIONS* const options,
					FILE* logstream)
{
	info(logstream, "checkout begin\n");
	
	/* check for NAN in any paramaters. */
	let _offset1 = advanfuncs->recordinfo.dataconfig->_offset1 ? 1 : 0;
	forcount(i, popmodel->ntheta) {
		let v = popmodel->theta[i];
		if (isnan(v)) {
			let n = i + _offset1;
			warning(logstream, "THETA(%i) is NAN\n", n);
		}
	}
	forcount(i, popmodel->nsigma) {
		let v = popmodel->sigma[i];
		if (isnan(v)) {
			let n = i + _offset1;
			warning(logstream, "SIGMA(%i) is NAN\n", n);
		}
	}
	forcount(i, popmodel->nomega) {
		forcount(j, i+1) {
			let v = popmodel->omega[i][j];
			if (isnan(v)) {
				let n1 = i + _offset1;
				let n2 = j + _offset1;
				warning(logstream, "OMEGA(%i,%i) is NAN\n", n1, n2);
			}
		}
	}

	SCATTEROPTIONS scatteroptions = { };
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

