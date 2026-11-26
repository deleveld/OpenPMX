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

#include "openpmx.h"
#include "idata.h"
#include "print.h"
#include "options.h"
#include "pmxstate.h"
#include "defines.h"
#include "advan/advan.h"
#include "utils/c22.h"
#include "utils/errctx.h"

#include <limits.h>

#include <gsl/gsl_cdf.h>

void pmx_npde(OPENPMX* source)
{
	info(0, "npde begin\n");

	/* ------------------------------------------------------------------ */
	/* private working copy – keeps the source ETAs / state untouched   */
	/* ------------------------------------------------------------------ */
	var pmx = pmx_copy(source);
	pmx_ensure_state(&pmx);

	/* ------------------------------------------------------------------ */
	/* simulations executed in working copy so that source ETA isnt modified
	 * Then the simulations are moved to source so they can be reused later */
	/* ------------------------------------------------------------------ */
	pmx_ensure_state(source);
	if (!source->state->idata.sim.dv) {
		SIMCONFIG simconfig = {
			.nsim = 1000
		};
		pmx_simulate(&pmx, &simconfig);

		/* move simulations to source, we dont need to free because we 
		 * already know there are none */ 
		source->state->idata.sim = pmx.state->idata.sim;
		memset(&pmx.state->idata.sim, 0, sizeof(pmx.state->idata.sim));
	} else
		info(0, "npde using previous %i simulations\n", source->state->idata.sim.nsim);

	let nsim  = source->state->idata.sim.nsim;
	let simdv = source->state->idata.sim.dv;

	/* allocate NPDE storage in the *source* so the caller can read it */
	var npde = idata_npde_alloc(&source->state->idata);

	/* scratch buffer for one observation’s simulated values */
	var buf = mallocvar(double, nsim);

	let advanfuncs = pmx.state->advanfuncs;
	let recordinfo = &advanfuncs->recordinfo;
	let dataconfig = recordinfo->dataconfig;
	let ndata      = recordinfo->ndata;
	var record     = dataconfig->records;
	let recordsize = dataconfig->recordfields.size;
	forcount(i, ndata) {
		if (RECORDINFO_EVID(recordinfo, record) != 0) {
			npde[i] = NAN;
			goto next;
		}

		let obs = RECORDINFO_DV(recordinfo, record);
		if (!isfinite(obs)) {
			npde[i] = NAN;
			goto next;
		}

		/* gather finite simulated values for this record */
		var m = 0;
		forcount(s, nsim) {
			let v = simdv[(size_t)s * ndata + i];
			if (isfinite(v))
				buf[m++] = v;
		}

		if (m == 0) {
			npde[i] = NAN;
		} else {
			/* mid-rank empirical probability */
			int below = 0, equal = 0;
			forcount(k, m) {
				if (buf[k] <  obs)
					below++;
				else if (buf[k] == obs)
					equal++;
			}
			var p = (below + 0.5 * equal) / m;

			/* clip to avoid ±Inf from the inverse normal */
			p = fmin(fmax(p, 0.5 / m), 1.0 - 0.5 / m);
			npde[i] = gsl_cdf_ugaussian_Pinv(p);
		}
next:
		record = RECORD_INDEX(record, recordsize, 1);
	}

	/* write output file */
	if (source->filename) {
		char fname[PATH_MAX];
		snprintf(fname, sizeof(fname), "%s" OPENPMX_NPDE, source->filename);
		var f = fopen(fname, "w");
		if (!f)
			fatal(0, "could not open npde file \"%s\"", fname);
		fprintf(f, OPENPMX_HEADER_FORMAT "\n", "NPDE");
		forcount(i, ndata) 
			fprintf(f, OPENPMX_TABLE_FORMAT "\n", npde[i]);
		fclose(f);
	}
	
	/* cleanup */
	free(buf);
	pmx_release_state(&pmx);

	info(0, "npde end\n");
}
