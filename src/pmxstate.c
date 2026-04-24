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
#include <string.h>

#include "openpmx.h"
#include "scatter.h"
#include "utils/c22.h"
#include "advan/advan.h"
#include "print.h"
#include "pmxstate.h"

static PMXSTATE* pmxstate_alloc(const OPENPMX* const pmx)
{
	let advanfuncs = advanfuncs_alloc(&pmx->data, &pmx->advan);
	ERRCTX errctx = { 0 };
	let popmodel = popmodel_init(pmx->theta, pmx->omega, pmx->sigma, &errctx);
	if (errctx.len)
		fatal(0, "%s: %s", __func__, errctx.errmsg);

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

void pmx_copy_popparam(OPENPMX* dest, const OPENPMX* const src)
{
	memcpy(dest->theta, src->theta, sizeof(dest->theta));
	memcpy(dest->sigma, src->sigma, sizeof(dest->sigma));
	memcpy(dest->omega, src->omega, sizeof(dest->omega));

	dest->result = (PMXRESULT) { 0 };
}

void pmx_cleanup(OPENPMX* pmx)
{
	pmxstate_free(pmx->state);
	pmx->state = 0;
}

void pmx_popmodel_writeback(OPENPMX* const pmx, const POPMODEL* const popmodel)
{
	let ntheta = popmodel->ntheta;
	forcount(i, ntheta) {
		pmx->theta[i].lower = popmodel->lower[i];
		pmx->theta[i].value = popmodel->theta[i];
		pmx->theta[i].upper = popmodel->upper[i];
		pmx->theta[i].type = popmodel->thetaestim[i];
	}

	let omega = popmodel->omega;
	let omegafixed = popmodel->omegafixed;
	var offset = 0;
	forcount(k, popmodel->nblock) {
		let ndim = popmodel->blockdim[k];
		let type = popmodel->blocktype[k];

		assert(pmx->omega[k].type == type);
		assert(pmx->omega[k].ndim == ndim);

		switch (type) {

			case OMEGA_DIAG: {
				forcount(i, ndim) {
					double v = omega[offset + i][offset + i];
					let f = omegafixed[offset + i][offset + i];
					/* fixed omegas on diagonal are negative */
					if (f != 0)
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
						/* fixed omegas on diagonal are negative */
						if (f != 0)
							v = -fabs(v);
						pmx->omega[k].values[blocki] = v;
						blocki++;
					}
				}
				break;
			}
			case OMEGA_SAME: {
				/* we dont need to do anything here */
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
		if (f != 0)
			v = -fabs(v);
		pmx->sigma[i] = v;
	}

	pmx->result = popmodel->result;
}

void pmx_set_theta(OPENPMX* dest, 
				   const int index,
				   typeof(((OPENPMX){0}).theta[0])* theta)
{
	let i = (dest->data._offset1) ? (index-1) : (index);
	assert(i >= 0);
	assert(i < OPENPMX_THETA_MAX);
	dest->theta[i] = *theta;
	
	dest->result = (PMXRESULT) { 0 };
}

