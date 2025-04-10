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
#include <assert.h>

#include "scatter.h"
#include "ievaluate.h"
#include "omegainfo.h"
#include "print.h"
#include "advan/advan.h"
#include "utils/c22.h"
#include "utils/various.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>

#include "openpmx_compile_options.h"

static void simulate_with_error_thread(INDIVID* const individ,
									   const ADVANFUNCS* const advanfuncs,
									   const POPMODEL* const popmodel,
									   const NONZERO* const nonzero,
									   const OPTIONS* const options,
									   const SCATTEROPTIONS* const scatteroptions)
{
	(void) nonzero;
	(void) options;
	
	let nomega = popmodel->nomega;
	double etacopy[OPENPMX_OMEGA_MAX] = { 0 };
	memcpy(etacopy, individ->eta, nomega * sizeof(double));

	let ievaluate_args = (IEVALUATE_ARGS) {
		.record = individ->record,
		.nrecord = individ->nrecord,
		.advanfuncs = advanfuncs,
		.popparam = popparam_init(popmodel, advanfuncs, etacopy),
		.logstream = scatteroptions->logstream,
	};

	struct timespec t3;
	clock_gettime(CLOCK_REALTIME, &t3);
	individual_simulate(&ievaluate_args,
						individ->imodel,
						individ->predictvars,
						individ->istate,
						individ->yhat,
						individ->yhatvar,
						individ->pred,
						individ->isimerr);
	individ->ineval += 1;
	individ->neval += 1;
	timespec_duration(&t3, &individ->eval_msec);
}

static void idata_resample_err(const IDATA* const idata,
							   const POPMODEL* const popmodel,
							   gsl_rng * const rng)
{
	let ndata = idata->ndata;
	let nsigma = popmodel->nsigma;

	/* fill in err values */
	var simerr = idata_alloc_simerr(idata);
	forcount(i, ndata) {
		let sigma = popmodel->sigma;
		forcount(j, nsigma) {
			let stddev = sqrt(sigma[j]);
			let v = gsl_ran_gaussian(rng, stddev);
			simerr[i * nsigma + j] = v;
		}
	}
}

static void idata_resample_eta(IDATA* const idata,
							   const POPMODEL* const popmodel,
							   gsl_rng * const rng)
{
	/* we need to know which rows and cols we need to resample */
	/* we cant do full matrix because some rows and coloums could be zero
	 * and we need to do resampling with cholesky whch would fail in that
	 * case due to non-positive-definite */
	/* this should fail for non-positive-definite because we could not sample
	 * from that matrix anyway and we should not sample from a different
	 * (i.e. corrected) one */
	var omegainfo = omegainfo_init(popmodel->nomega, popmodel->omega, popmodel->omegafixed);

	/* first zero everything */
	let nomega = idata->nomega;
	forcount(i, idata->nindivid) {
		let individ = &idata->individ[i];

		individ->obs_min2ll = 0.;
		individ->obs_lndet = 0.;
		individ->eta_min2ll = 0.;
		individ->icov_lndet = 0.;
		individ->iobjfn = 0.;

		let icov = individ->icov;
		memset(icov, 0, nomega * nomega * sizeof(double));

		let eta = individ->eta;
		memset(eta, 0, nomega * sizeof(double));
	}

	/* maybe we have no etas to sample */
	let n = omegainfo.nonzero.n;
	if (n == 0)
		return;

	/* we have to make cholesky here for resampling */
	var cholesky = gsl_matrix_view_array(omegainfo.nonzero.choleskydata, n, n);

	/* resample eta matrix */
	var v = gsl_vector_alloc(n);
	var s = gsl_vector_alloc(n);
	forcount(i, idata->nindivid) {
		let individ = &idata->individ[i];

		/* multiply random N(0,1) vector by cholesky matrix */
		forcount(j, v->size) {
			let x = gsl_ran_gaussian(rng, 1.);
			gsl_vector_set(v, j, x);
		}
		gsl_blas_dgemv(CblasNoTrans, 1., &cholesky.matrix, v, 0., s);

		/* fill in the non-zero places in eta */
		let eta = individ->eta;
		let nonzero = omegainfo.nonzero.rowcol;
		forcount(j, omegainfo.nonzero.n) {
			let c = nonzero[j];
			let x = gsl_vector_get(s, j);
			eta[c] = x;
		}
	}
	gsl_vector_free(v);
	gsl_vector_free(s);
}

/* after idata_predict_dv then DV in thet dataset is replaced by observation
 * simulation with noise. pred is set to zero and yhat is the observation
 * without noise. If the error has not been resampled yet by calling
 * idata_resample_err then the error will be zero*/
static void idata_predict_dv(IDATA* const idata,
							const ADVANFUNCS* const advanfuncs,
							const POPMODEL* popmodel,
							const OPTIONS* const options)
{
	/* if error has not been resampled yet, then use zero for error */
	idata_alloc_simerr(idata);

	/* simulation writes into pred (with error) and yhat (without error) */
	/* TODO: we have to not write non-DV objects in simulation. You mean predict?
	 * Im not sure I understand anymore */
	SCATTEROPTIONS scatteroptions = { 0 };
	scatter_threads(idata, advanfuncs, popmodel, 0, options, &scatteroptions, simulate_with_error_thread);

	/* make sure data is writable */
	var recordinfo = &advanfuncs->recordinfo;
	var dataconfig = recordinfo->dataconfig;
	var writeable = dataconfig->writeable;
	var readabledata = dataconfig->records;
	if (writeable == 0)
		fatal(0, "data not writeable\n");
	if (writeable != readabledata)
		fatal(0, "writeable data does not match readable data\n");

	var offsetDV = recordinfo->offsetDV;
	var ptrdv = (double*)((char*)writeable + offsetDV);
	let pred = idata->individ[0].pred;

	/* write back from pred (obs *with* noise) into DV
	 * and zero pred after this */
	let ndata = recordinfo->ndata;
	var recordsize = dataconfig->recordfields.size;
	forcount(i, ndata) {
		*ptrdv = pred[i];
		ptrdv = (double*)((char*)ptrdv + recordsize);
		pred[i] = 0.;
	}
}

static void idata_resample(IDATA* const idata,
						   const POPMODEL* popmodel,
						   gsl_rng* rng)
{
	idata_resample_err(idata, popmodel, rng);
	idata_resample_eta(idata, popmodel, rng);
}

void pmx_simulate(OPENPMX* pmx, const SIMCONFIG* const simconfig)
{
	pmxstate_ensure(pmx);
	var pstate = pmx->state;

	var popmodel = popmodel_init(pmx);
	var options = options_init(pmx);
	if (simconfig)
		options.simulate = simconfig_default(simconfig);

	/* allocate random number generator */
	if (!pstate->rng) {
		var seed = options.simulate.seed;
		pstate->rng = gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(pstate->rng, seed);
	}
	assert(pstate->rng);

	idata_resample(&pstate->idata, &popmodel, pstate->rng);
	idata_predict_dv(&pstate->idata, pstate->advanfuncs, &popmodel, &options);
}

