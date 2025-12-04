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

#include <assert.h>
#include <float.h>
#include <unistd.h>

#include "openpmx.h"
#include "scatter.h"
#include "utils/c22.h"
#include "advan/advan.h"
#include "print.h"
#include "pmxstate.h"

static PMXSTATE* pmxstate_alloc(const OPENPMX* const pmx)
{
	let advanfuncs = advanfuncs_alloc(&pmx->data, &pmx->advan);
	let popmodel = popmodel_init(pmx);

	var temp = (PMXSTATE) {
		.advanfuncs = advanfuncs,
		.idata = idata_construct(&advanfuncs->recordinfo,
								 popmodel.ntheta,
								 popmodel.nomega,
								 popmodel.nsigma,
								 advanfuncs->nstate,
								 advanfuncs->advanconfig->imodelfields.size,
								 advanfuncs->advanconfig->predictfields.size),
		.tablecount = 0,
		.rng = 0,
	};

	/* have to allocate because its a opaque pointer so we have to memcpy */
	var ret = mallocvar(PMXSTATE, 1);
	memcpy(ret, &temp, sizeof(temp));
	return ret;
}

void pmxstate_ensure(OPENPMX* const pmx)
{
	if (!pmx->state)
		pmx->state = pmxstate_alloc(pmx);
}

static void pmxstate_free(PMXSTATE* pstate)
{
	scatter_cleanup();

	advanfuncs_free((void*)pstate->advanfuncs); /* keeping advanfuncs const in PMXSTATE is nice, but this cludge is needed to free */
	idata_destruct(&pstate->idata);

	if (pstate->rng)
		gsl_rng_free(pstate->rng);

	free(pstate);
}

OPENPMX pmx_copy(const OPENPMX* const pmx)
{
	OPENPMX ret = *pmx;
	ret.state = 0;
	return ret;
}

void pmx_copy_popparams(OPENPMX* dest, const OPENPMX* const src)
{
	memcpy(dest->theta, src->theta, sizeof(dest->theta));
	memcpy(dest->sigma, src->sigma, sizeof(dest->sigma));
	memcpy(dest->omega, src->omega, sizeof(dest->omega));

	dest->result = (PMXRESULT) { .objfn = DBL_MAX,
								 .type = OBJFN_INVALID,
								 .nparam = 0,
								 .neval = 0 };
}

void pmx_cleanup(OPENPMX* pmx)
{
	pmxstate_free(pmx->state);
	pmx->state = 0;
}

void pmx_update_from_popmodel(OPENPMX* const pmx, const POPMODEL* const popmodel)
{
	let theta = popmodel->theta;
	let ntheta = popmodel->ntheta;
	let thetaestim = popmodel->thetaestim;
	forcount(i, ntheta) {
		if (thetaestim[i] == ESTIMATE)
			pmx->theta[i].value = theta[i];
	}

	let omega = popmodel->omega;
	let omegafixed = popmodel->omegafixed;
	var offset = 0;
	forcount(k, popmodel->nblock) {
		let ndim = popmodel->blockdim[k];
		let type = popmodel->blocktype[k];

		assert((int)pmx->omega[k].type == type);
		assert(pmx->omega[k].ndim == ndim);

		switch (type) {

			case OMEGA_DIAG: {
				forcount(i, ndim) {
					double v = omega[offset + i][offset + i];
					let f = omegafixed[offset + i][offset + i];
					if (f != 0 && v != 0.)
						v = -fabs(v);
					pmx->omega[k].values[i] = v;
				}
				break;
			}
			case OMEGA_BLOCK: {
				int blocki = 0;
				forcount(i, ndim) {
					forcount(j, i + 1) {
						double v = omega[offset + i][offset + j];
						let f = omegafixed[offset + i][offset + j];
						if (f != 0 && v != 0.)
							v = -fabs(v);
						pmx->omega[k].values[blocki] = v;
						blocki++;
					}
				}
				break;
			}
			case OMEGA_SAME: {
				/* we dont need to do anything here since the block is already */
				break;
			}
			default:
				fatal(0, "Invalid block type (%i) if size %i\n", type, ndim);
				break;
		}
		offset += ndim;
	}

	let sigma = popmodel->sigma;
	let nsigma = popmodel->nsigma;
	let sigmafixed = popmodel->sigmafixed;
	forcount(i, nsigma) {
		var v = sigma[i];
		var f = sigmafixed[i];
		if (f == 0)
			pmx->sigma[i] = v;
	}

	pmx->result = popmodel->result;
}

