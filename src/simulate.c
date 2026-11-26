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

/// This file implements the top-level function for simulation.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "simulate.h"
#include "scatter.h"
#include "ievaluate.h"
#include "omegainfo.h"
#include "print.h"
#include "pmxstate.h"
#include "advan/advan.h"
#include "utils/c22.h"
#include "utils/various.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>

void simulate_with_error_thread(INDIVID* const individ,
								const ADVANFUNCS* const advanfuncs,
								const POPMODEL* const popmodel,
								const NONZERO* const nonzero,
								const OPTIONS* const options,
								const SCATTEROPTIONS* const scatteroptions)
{
	(void) nonzero;
	(void) options;
	
	let nomega = popmodel->nomega;
	double etacopy[OPENPMX_OMEGA_MAX];
	forcount(i, OPENPMX_OMEGA_MAX)
		etacopy[i] = NAN;
	memcpy(etacopy, individ->eta, nomega * sizeof(double));

	struct timespec t3;
	clock_gettime(CLOCK_MONOTONIC, &t3);
	let ievaluate_args = ievaluate_args_init(individ->record,
											 individ->nrecord,
											 advanfuncs,
											 popmodel->theta,
											 popmodel->ntheta,
											 etacopy,
											 popmodel->nomega,
											 popmodel->sigma,
											 popmodel->nsigma,
											 scatteroptions->logstream);
	individual_simulate(&ievaluate_args,
						individ->imodel,
						individ->predictvars,
						individ->istate,
						individ->yhat,
						individ->yhatvar,
						individ->pred,
						individ->simerr);
	individ->ineval += 1;
	timespec_duration(&t3, &individ->eval_msec);
}

/// For simulation random observation residuals are chosen according to
/// sigma. 
void idata_resample_err(const IDATA* const idata,
						const ADVANFUNCS* const advanfuncs,
						const POPMODEL* const popmodel,
						gsl_rng * const rng)
{
	let ndata = idata->ndata;
	let nsigma = popmodel->nsigma;

	var recordinfo = &advanfuncs->recordinfo;
	let dataconfig = recordinfo->dataconfig;
	let recordsize = dataconfig->recordfields.size;
	let sigma = popmodel->sigma;

	/* fill in err values, only observations get random error */
	var simerr = idata_simerr_alloc(idata);
	const RECORD* record = dataconfig->records;
	forcount(i, ndata) {
		if (RECORDINFO_EVID(recordinfo, record) == 0) {
			forcount(j, nsigma) {
				let stddev = sqrt(sigma[j]);
				simerr[i * nsigma + j] = gsl_ran_gaussian(rng, stddev);
			}
		} else {
			forcount(j, nsigma)
				simerr[i * nsigma + j] = 0.;
		}
		record = (const RECORD*)((const char*)record + recordsize);
	}
}

/// For simulation random ETA values are chosen according to the omega
/// matrix.
void idata_resample_eta(IDATA* const idata,
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

/// Upon simulation the individualized values the objective function components are set to zero.
	let nomega = idata->nomega;
	forcount(i, idata->nindivid) {
		let individ = &idata->individ[i];

		individ->obs_min2ll = 0.;
		individ->obs_lndet = 0.;
		individ->obs_logp = 0.;
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

/// Upon simulation the individualized eta values are randomly selected from 
/// from the distribution of the omega matrix.
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
//		gsl_blas_dgemv(CblasNoTrans, 1., &cholesky.matrix, v, 0., s);
		gsl_blas_dcopy(v, s);  
		gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, &cholesky.matrix, s);

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

/// After simulation a prediction occurs and the DV in the dataset is
/// replaced by observation with noise. PRED is set to zero and yhat is
/// the observation without noise. 
static void idata_predict_dv(IDATA* const idata,
							 const ADVANFUNCS* const advanfuncs,
							 const POPMODEL* popmodel,
							 const OPTIONS* const options)
{
	/* make sure error has been allocated. */
	assert(idata->individ[0].simerr);

	/* simulation writes into pred (with error) and yhat (without error) */
	SCATTEROPTIONS scatteroptions = { };
	scatter_threads(idata, advanfuncs, popmodel, 0, options, &scatteroptions, simulate_with_error_thread);
}

static void idata_writeback_dv(IDATA* const idata,
							   const ADVANFUNCS* const advanfuncs)
{
	/* make sure data is writable */
	var recordinfo = &advanfuncs->recordinfo;
	var dataconfig = recordinfo->dataconfig;
	var writeable = dataconfig->writeable;
	var readabledata = dataconfig->records;
	if (writeable == 0)
		fatal(0, "data not writeable\n");
	if (writeable != readabledata)
		fatal(0, "writeable data does not match readable data\n");

	let offsetDV = recordinfo->offsetDV;
	let pred = idata->individ[0].pred;
	assert(offsetDV >= 0);

	/* write back from pred (obs *with* noise) into DV
	 * and zero pred after this */
	let ndata = recordinfo->ndata;
	let recordsize = dataconfig->recordfields.size;
	RECORD* record = writeable;
	forcount(i, ndata) {
		if (RECORDINFO_EVID(recordinfo, record) == 0)
			*(double*)((char*)record + offsetDV) = pred[i];
		pred[i] = 0.;
		record = (RECORD*)((char*)record + recordsize);
	}
}

void pmx_simulate(OPENPMX* pmx, const SIMCONFIG* const simconfig)
{
	pmx_ensure_state(pmx);
	var pstate = pmx->state;

	var options = options_init(pmx);
	var nsim = 0;
	if (simconfig) {
		options.simulate = simconfig_default(simconfig);
		nsim = simconfig->nsim;
	}

	if (nsim < 0) 
		fatal(0, "sim output nsim < 0\n");

	ERRCTX errctx = { 0 };
	var popmodel = popmodel_init(pmx->theta, pmx->omega, pmx->sigma, &errctx);
	if (errctx.len)
		fatal(0, "%s", errctx.errmsg);
	
	pmx_ensure_state_rng(pmx, &options);

	/// If nsim == 0 then only one simulation is done and the results
	/// written back to the dataset DV.
	var idata = &pstate->idata;
	var advanfuncs = pstate->advanfuncs;
	if (nsim == 0) {
		info(0, "simulate begin modify DV\n");

		idata_sim_free(idata);
		idata_resample_eta(idata, &popmodel, pstate->rng);
		idata_resample_err(idata, advanfuncs, &popmodel, pstate->rng);
		idata_predict_dv(idata, advanfuncs, &popmodel, &options);
		
		/* write-back into dataset */
		idata_writeback_dv(idata, advanfuncs);

	/// If nsim > 0 then only many copies of the simulation are done and
	/// the result DV values saved in the individual data.
	} else {
		var t = 0.;
		var next_time = 0.;
		struct timespec begin;
		clock_gettime(CLOCK_MONOTONIC, &begin);

		info(0, "sim begin %i\n", nsim);

		var recordinfo = &advanfuncs->recordinfo;
		let ndata = recordinfo->ndata;
		var dv = idata_sim_alloc(idata, nsim);
		forcount(i, nsim) {
			struct timespec now;
			clock_gettime(CLOCK_MONOTONIC, &now);
			t = timespec_time_difference(&begin, &now) / 1000.;
			if (i == 0 || t >= next_time) {
				info(0, "time %.3f simulate %i\n", next_time, i+1);
				next_time = t + 1.;
			} 
			idata_resample_eta(idata, &popmodel, pstate->rng);
			idata_resample_err(idata, advanfuncs, &popmodel, pstate->rng);
			idata_predict_dv(idata, advanfuncs, &popmodel, &options);

			/* collect DV in pointer and zero pred */
			let pred = idata->individ[0].pred;
			memcpyvar(&dv[i * ndata], pred, ndata);
			memset(pred, 0, ndata*sizeof(double));
		}
		info(0, "time %.3f simulate %i\n", t, nsim);
	}
	info(0, "sim end\n");
}
