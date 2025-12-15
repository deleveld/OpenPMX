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
#include <float.h>

#include "idata.h"
#include "defines.h"
#include "utils/c22.h"
#include "utils/various.h"

#include <gsl/gsl_math.h>

IDATA idata_construct(const RECORDINFO* const recordinfo,
					  const int ntheta,
					  const int nomega,
					  const int nsigma,
					  const int nstate,
					  const int imodel_size,
					  const int predictvars_size)
{
	let data = recordinfo->dataconfig->records;
	let ndata = recordinfo->ndata;
	let nindivid = recordinfo->nindivid;
	let nobs = recordinfo->nobs;

	assert(ntheta <= OPENPMX_THETA_MAX);
	assert(nomega <= OPENPMX_OMEGA_MAX);
	assert(nsigma <= OPENPMX_SIGMA_MAX);
	assert(nstate <= OPENPMX_STATE_MAX);

	var imodel = (IMODEL*)calloc(ndata, imodel_size);
	var state = callocvar(double, ndata * nstate);
	var eta = callocvar(double, nindivid * nomega);
	var icov = callocvar(double, nindivid * nomega * nomega);
	var yhat = callocvar(double, ndata);
	var yhatvar = callocvar(double, ndata);
	var pred = callocvar(double, ndata);
	var predictvars = (PREDICTVARS*)calloc(ndata, predictvars_size);

	/* fill in individual information */
	var individ = mallocvar(INDIVID, nindivid);
	var n = 0;
	var i = 0;
	while (i < ndata) {
		var nobsi = 0;
		var nrecord = 0;
		let thisid = RECORDINFO_ID(recordinfo, RECORDINFO_INDEX(recordinfo, data, i));
		while (i+nrecord < ndata && RECORDINFO_ID(recordinfo, RECORDINFO_INDEX(recordinfo, data, i+nrecord)) == thisid) {
			if (RECORDINFO_EVID(recordinfo, RECORDINFO_INDEX(recordinfo, data, i+nrecord)) == 0)
				++nobsi;
			++nrecord;
		}

		/* initialise INDIVID for each individual. We have to do this in place
		 * and then copy over for the const values, otherwise we are copy over a
		 * read-only location. This is explicitly allowed in C to write over
		 * const struct members if the memory has been malloced. */
		var temp = (INDIVID) {
			.ID = thisid,
			.record = RECORDINFO_INDEX(recordinfo, data, i),
			.nrecord = nrecord,
			.nobs = nobsi,

			/* individual data are views into the main matrix */
			.imodel = (IMODEL*)((char*)imodel + i * imodel_size),
			.istate = &state[i * nstate],
			.eta = &eta[n * nomega],
			.icov = &icov[n * nomega * nomega],
			.yhat = &yhat[i],
			.yhatvar = &yhatvar[i],
			.pred = &pred[i],
			.predictvars = (PREDICTVARS*)((char*)predictvars + i * predictvars_size),
			.icovsample = 0,
			.icovweight = 0,
			.isimerr = 0,

			.obs_min2ll = 0.,
			.obs_lndet = 0.,
			.eta_min2ll = 0.,
			.icov_lndet = 0.,
			.iobjfn = DBL_MAX,

			.eval_msec = 1. + nrecord,		/* first guess as to evaluation time is the number of records */
			.stage1_msec = 1. + nrecord,
			.ineval = 0,
		};
		var iptr = &individ[n];
		memcpy(iptr, &temp, sizeof(temp));

		++n;
		i += nrecord;
	}

	return (IDATA) {
		.ndata = ndata,
		.nobs = nobs,
		.nindivid = nindivid,
		.ntheta = ntheta,
		.nomega = nomega,
		.nsigma = nsigma,
		.nstate = nstate,
		.imodel_size = imodel_size,
		.predictvars_size = predictvars_size,

		.individ = individ,
	};
}

void idata_destruct(IDATA* const idata)
{
	/* first individual has the all the memory */
	var firstindivid = &idata->individ[0];
	free(firstindivid->imodel);
	free(firstindivid->istate);
	free(firstindivid->eta);
	free(firstindivid->icov);
	free(firstindivid->yhat);
	free(firstindivid->yhatvar);
	free(firstindivid->pred);
	free(firstindivid->predictvars);
	free(firstindivid->icovsample);	/* could be zero. allocated in idata_alloc_icovridata */
	free(firstindivid->icovweight);	/* could be zero. allocated in idata_alloc_icovridata */
	free(firstindivid->isimerr);	/* could be zero. allocated in idata_alloc_simerr */

	/* now we can delete the individ array */
	free(idata->individ);
}

double* idata_alloc_simerr(const IDATA* const idata)
{
	let individ = idata->individ;
	if (individ[0].isimerr == 0) {
		let ndata = idata->ndata;
		let nsigma = idata->nsigma;
		var simerr = mallocvar(double, ndata * nsigma);

		let nindivid = idata->nindivid;
		var ioffset = 0;
		forcount(i, nindivid) {
			individ[i].isimerr = &simerr[ioffset * nsigma];
			ioffset += individ[i].nrecord;
		}
	}
	return individ[0].isimerr;
}

void idata_free_simerr(const IDATA* const idata)
{
	let individ = idata->individ;
	let nindivid = idata->nindivid;
	var simerr = individ[0].isimerr;
	forcount(i, nindivid)
		individ[i].isimerr = 0;
	free(simerr);
}

void idata_alloc_icovresample(const IDATA* const idata)
{
	let individ = idata->individ;
	let nindivid = idata->nindivid;
	let nomega = idata->nomega;
	if (individ[0].icovsample == 0) {
		assert(individ[0].icovweight == 0);
		var icovsample = callocvar(double, nindivid * 2 * nomega * nomega);
		var icovweight = callocvar(double, nindivid * 2 * nomega);

		forcount(i, nindivid) {
			individ[i].icovsample = &icovsample[i * 2 * nomega * nomega];
			individ[i].icovweight = &icovweight[i * 2 * nomega];
		}
	} else {
		forcount(i, nindivid * 2 * nomega * nomega)
			individ[0].icovsample[i] = 0.;
		forcount(i, nindivid * 2 * nomega)
			individ[0].icovweight[i] = 0.;
	}
}

void idata_free_icovresample(const IDATA* const idata)
{
	let individ = idata->individ;
	let nindivid = idata->nindivid;
	var icovsample = individ[0].icovsample;
	var icovweight = individ[0].icovweight;
	forcount(i, nindivid) {
		individ[i].icovsample = 0;
		individ[i].icovweight = 0;
	}
	free(icovsample);
	free(icovweight);
}

int idata_ineval(const IDATA* const idata, const bool reset)
{
	let individ = idata->individ;
	let nindivid = idata->nindivid;
	var ineval = 0;
	forcount(i, nindivid) {
		ineval = individ[i].ineval;
		if (reset)
			individ[i].ineval = 0;
	}
	return ineval;
}

void idata_reset_eta(IDATA* const idata, const double* eta)
{
	assert(eta);
	var firstindivid = &idata->individ[0];
	memcpy(firstindivid->eta, eta, idata->nindivid * idata->nomega * sizeof(double));
}

void table_phi_idata(const char* filename,
					 const IDATA* const idata,
					 const bool _offset1)
{
	assert(filename);
	var f = results_fopen(filename, OPENPMX_PHIFILE, "w");
	assert(f);

	let nomega = idata->nomega;
	let indexoffset = _offset1 ? 1 : 0;

	fprintf(f, OPENPMX_SFORMAT, "SUBJECT_NO");
	fprintf(f, OPENPMX_HEADER_FORMAT, "ID");
	forcount(i, nomega) {
		char temp[128];
		sprintf(temp, "ETA(%i)", i + indexoffset);
		fprintf(f, OPENPMX_HEADER_FORMAT, temp);
	}
	forcount(i, nomega) {
		forcount(j, i+1) {
			char temp[128];
			sprintf(temp, "ETC(%i,%i)", i + indexoffset, j + indexoffset);
			fprintf(f, OPENPMX_HEADER_FORMAT, temp);
		}
	}
	fprintf(f, OPENPMX_HEADER_FORMAT, "OBJ");
	
	fprintf(f, OPENPMX_HEADER_FORMAT, "TEVAL");
	fprintf(f, OPENPMX_HEADER_FORMAT, "TSTAGE1");
	fprintf(f, OPENPMX_SFORMAT, "INEVAL");
	fprintf(f, "\n");

	forcount(k, idata->nindivid) {
		fprintf(f, OPENPMX_IFORMAT, k + 1);

		let individ = &idata->individ[k];
		let id = individ->ID;
		fprintf(f, OPENPMX_TABLE_FORMAT, id);

		var eta = individ->eta;
		forcount(i, nomega) {
			let v = eta[i];
			fprintf(f, OPENPMX_TABLE_FORMAT, v);
		}
		let icov = individ->icov;
		forcount(i, nomega) {
			forcount(j, i+1) {
				let v = icov[i * nomega + j];
				fprintf(f, OPENPMX_TABLE_FORMAT, v);
			}
		}

		let iobjfn = individ->iobjfn;
		fprintf(f, OPENPMX_TABLE_FORMAT, iobjfn);

		let teval = individ->eval_msec;
		let tstage1 = individ->stage1_msec;
		let ineval = individ->ineval;
		fprintf(f, OPENPMX_TABLE_FORMAT, teval);
		fprintf(f, OPENPMX_TABLE_FORMAT, tstage1);
		fprintf(f, OPENPMX_IFORMAT, ineval);

		fprintf(f, "\n");
	}

	fclose(f);
}

void table_icov_resample_idata(const char* filename,
							   const IDATA* const idata,
							   const bool _offset1)
{
	assert(filename);
	assert(idata->individ[0].icovweight != 0);

	var f = results_fopen(filename, OPENPMX_ICOVRESAMPLEFILE, "w");
	assert(f);

	let nomega = idata->nomega;
	let indexoffset = _offset1 ? 1 : 0;

	fprintf(f, OPENPMX_SFORMAT, "SUBJECT_NO");
	fprintf(f, OPENPMX_HEADER_FORMAT, "ID");
	fprintf(f, OPENPMX_HEADER_FORMAT, "RESAMPLE");
	fprintf(f, OPENPMX_HEADER_FORMAT, "WEIGHT");
	forcount(i, nomega) {
		char temp[128];
		sprintf(temp, "ETA(%i)", i + indexoffset);
		fprintf(f, OPENPMX_HEADER_FORMAT, temp);
	}
	fprintf(f, "\n");

	forcount(k, idata->nindivid) {
		let individ = &idata->individ[k];
		let id = individ->ID;
		let icovweight = individ->icovweight;
		let icovsample = individ->icovsample;

		let eta = individ->eta;
		fprintf(f, OPENPMX_IFORMAT, k + 1);
		fprintf(f, OPENPMX_TABLE_FORMAT, id);
		fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
		fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
		forcount(j, nomega) {
			let v = eta[j];
			fprintf(f, OPENPMX_TABLE_FORMAT, v);
		}
		fprintf(f, "\n");

		forcount(i, 2 * nomega) {
			let wgt = icovweight[i];
			if (wgt > 0.) {
				fprintf(f, OPENPMX_IFORMAT, k + 1);
				fprintf(f, OPENPMX_TABLE_FORMAT, id);
				fprintf(f, OPENPMX_TABLE_FORMAT, i + 1.);
				fprintf(f, OPENPMX_TABLE_FORMAT, wgt);
				forcount(j, nomega) {
					let v = icovsample[i * nomega + j];
					fprintf(f, OPENPMX_TABLE_FORMAT, v);
				}
				fprintf(f, "\n");
			}
		}
	}
	fclose(f);
}

double idata_objfn(const IDATA* const idata,
				   const double omega_nonzero_lndet)
{
	let nindivid = idata->nindivid;

	/* doing the sum by type makes sure we are adding numbers of comparable
	 * magnitude which helps with accuracy */
	double objfn1 = 0.;
	double objfn2 = 0.;
	double objfn3 = 0.;
	double objfn4 = 0.;
	double objfn5 = 0.;
	forcount(k, nindivid) {
		let individ = &idata->individ[k];

		let term1 = individ->obs_lndet;
		let term2 = individ->obs_min2ll;
		let term3 = individ->eta_min2ll;
		var term4 = omega_nonzero_lndet;
 		let term5 = individ->icov_lndet;
		if (individ->nobs == 0) {
			assert(term1 == 0.);
			assert(term2 == 0.);
			assert(term3 == 0.);
			term4 = 0.;				/* population term does not count if no observations in the individual */
			assert(term5 == 0.);
		}

		assert(gsl_finite(term1) == 1);
		assert(gsl_finite(term2) == 1);
		assert(gsl_finite(term3) == 1);
		assert(gsl_finite(term4) == 1);
		assert(gsl_finite(term5) == 1);

		let iobjfn = term1 + term2 + term3 + term4 + term5;
		individ->iobjfn = iobjfn;

		objfn1 += term1;
		objfn2 += term2;
		objfn3 += term3;
		objfn4 += term4;
		objfn5 += term5;
	}
	let objfn = objfn1 + objfn2 + objfn3 + objfn4 + objfn5;
	assert(gsl_finite(objfn) == 1);
	return objfn;
}



