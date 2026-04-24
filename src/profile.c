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

#include <limits.h>

#include "openpmx.h"
#include "popmodel.h"
#include "pmxstate.h"
#include "omegafixed.h"
#include "print.h"
#include "utils/c22.h"
#include "utils/errctx.h"

/// This file implements a function that allows likelihood profiling.
/// It is available in openpmxtran as `profile()`.

OPENPMX pmx_profile_evaluate(const OPENPMX* const source, PROFILECONFIG* const args)
{
	ERRCTX errctx = { 0 };

	if (args->type == PROFILE_INVALID) 
		fatal(0, "%s: profile type invalid\n", __func__);

	/* make a POPMODEL and set the profile paramaters */
	var popmodel = popmodel_init(source->theta, source->omega, source->sigma, &errctx);
	if (errctx.len)
		fatal(0, "%s: %s", __func__, errctx.errmsg);
		
	const char* profile_type = "unknown";
	var index = args->index - (source->data._offset1 ? 1 : 0);
	if (args->type == PROFILE_THETA) {
		profile_type = "theta";
		if (index < 0 || index >= popmodel.ntheta)
			fatal(0, "%s: profile %s index (%i) out of bounds", __func__, profile_type, args->index);
		popmodel.upper[index] = args->value;
		popmodel.theta[index] = args->value;
		popmodel.lower[index] = args->value;
		popmodel.thetaestim[index] = FIXED;

	} else if (args->type == PROFILE_OMEGA) {
		profile_type = "omega";
		if (index < 0 || index >= popmodel.nomega)
			fatal(0, "%s: profile %s index (%i) out of bounds", __func__, profile_type, args->index);
		if (popmodel.omegafixed[index][index] == OMEGAFIXED_SAME)
			fatal(0, "%s: cannot profile part of omega same block\n", __func__);
		popmodel.omega[index][index] = args->value;
		popmodel.omegafixed[index][index] = OMEGAFIXED_FIXED;

	} else if (args->type == PROFILE_SIGMA) {
		profile_type = "sigma";
		if (index < 0 || index >= popmodel.nsigma)
			fatal(0, "%s: profile %s index (%i) out of bounds", __func__, profile_type, args->index);
		popmodel.sigma[index] = args->value;
		popmodel.sigmafixed[index] = 1;

	} else 
		fatal(0, "%s: profile type (%i) invalid", __func__, args->type);

	/* make our OPENPMX object by taking the source and reset state and
	 * copy the new settings to it */
	/* make a filename to write results to */
	var ret = *source;
	ret.state = 0;
	char filename[PATH_MAX] = { 0 };
	if (source->filename) {
		snprintf(filename, sizeof(filename), "%s.profile.evaluate.%s.%i", source->filename, profile_type, args->index);
		ret.filename = filename;
	}
	pmx_popmodel_writeback(&ret, &popmodel);
	
	/* do the estimation or evaluation */
	if (args->estimconfig.maxeval == 1)
		pmx_evaluate(&ret, &args->estimconfig.stage1);
	else 
		pmx_estimate(&ret, &args->estimconfig);
		
	/* dont let the filename dangle */
	ret.filename = 0;
	
	return ret;
}
