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
#include <math.h>
#include <float.h>

#include "popmodel.h"
#include "print.h"
#include "defines.h"
#include "utils/c22.h"
#include "utils/various.h"

POPMODEL popmodel_init(const OPENPMX* const pmx)
{
	let theta = pmx->theta;
	let omegablocks = pmx->omega;
	let sigma = pmx->sigma;
	
	POPMODEL ret = { };
	ret.ntheta = 0;
	ret.nblock = 0;
	ret.nomega = 0;
	ret.nsigma = 0;

	/* initialize to invalid values to catch errors */
	forcount(i, OPENPMX_THETA_MAX) {
		ret.lower[i] = NAN;
		ret.theta[i] = NAN;
		ret.upper[i] = NAN;
		ret.thetaestim[i] = THETA_INVALID;
	}
	forcount(i, OPENPMX_OMEGABLOCK_MAX) {
		ret.blocktype[i] = OMEGA_INVALID;
		ret.blockdim[i] = 0;
	}
	forcount(i, OPENPMX_OMEGA_MAX) 
		forcount(j, OPENPMX_OMEGA_MAX)
			ret.omega[i][j] = NAN;
	forcount(i, OPENPMX_SIGMA_MAX) {
		ret.sigma[i] = NAN;
		ret.sigmafixed[i] = 1;
	}

	/* count thetas, omegas, sigmas */
	forcount(i, OPENPMX_THETA_MAX) {
		let est = theta[i].type;
		if (est == THETA_INVALID) 
			break;
		ret.ntheta = i + 1;

		if (!isfinite(theta[i].lower) ||
			!isfinite(theta[i].value) ||
			!isfinite(theta[i].upper))
			fatal(0, "THETA invalid { %g, %g, %g }\n", theta[i].lower, theta[i].value, theta[i].upper);
		
		ret.lower[i] = theta[i].lower;
		ret.theta[i] = theta[i].value;
		ret.upper[i] = theta[i].upper;

		if (theta[i].value < theta[i].lower ||
			theta[i].value > theta[i].upper)
			fatal(0, "THETA outside boundary { %g, %g, %g }\n", theta[i].lower, theta[i].value, theta[i].upper);

		if (est == ESTIMATE) {
			if (theta[i].value == theta[i].lower ||
				theta[i].value == theta[i].upper)
				fatal(0, "THETA at boundary { %g, %g, %g, ESTIMATE }\n", theta[i].lower, theta[i].value, theta[i].upper);
		}
		ret.thetaestim[i] = est;
	}
	for(int i=ret.ntheta; i<OPENPMX_THETA_MAX; i++) {
		if (theta[i].lower != 0. ||
			theta[i].value != 0. ||
			theta[i].upper != 0. ||
			theta[i].type != THETA_INVALID)
			fatal(0, "THETA not fully initialized { %g, %g, %g, ??? }\n", theta[i].lower, theta[i].value, theta[i].upper);
	}

	/* find the dimensions to the whole omega matrix by adding the size of each block */
	forcount(i, OPENPMX_OMEGABLOCK_MAX) {
		let type = omegablocks[i].type;
		if (type == OMEGA_INVALID)
			break;
		
		let ndim = omegablocks[i].ndim;
		ret.blocktype[i] = type;
		ret.blockdim[i] = ndim;
		ret.nomega += ndim;
		ret.nblock += 1;
	}
	if (ret.nomega > OPENPMX_OMEGA_MAX)
		fatal(0, "Omega size (%i) is too large, max is %i\n", ret.nomega, OPENPMX_OMEGA_MAX);

	/* initialize the omega matrix from the blocks */
	forcount(i, ret.nomega) 
		forcount(j, ret.nomega)
			ret.omega[i][j] = 0.;
	var d = 0;
	forcount(i, ret.nblock) {
		let ndim = omegablocks[i].ndim;
		let v = omegablocks[i].values;
		let type = omegablocks[i].type;
		var n = 0;
		if (type == OMEGA_BLOCK) {
			forcount(r, ndim) {
				/* off diagonal */
				for (var c=0; c<r; c++) {
					let value = v[n];
					ret.omega[d + c][d + r] = value;
					ret.omega[d + r][d + c] = value;
					++n;
				}
				/* on diagonal */
				ret.omega[d + r][d + r] = v[n];
				++n;
			}

		/* get entries from elsewhere set omega elements */
		} else if (type == OMEGA_SAME) {
			forcount(r, ndim) {
				for (var c=0; c<r; c++) {
					let value = ret.omega[d + c - ndim][d + r - ndim];
					ret.omega[d + c][d + r] = value;
					ret.omega[d + r][d + c] = value;

					ret.omegafixed[d + c][d + r] = 2;
					ret.omegafixed[d + r][d + c] = 2;
				}
				let value = ret.omega[d + r - ndim][d + r - ndim];
				ret.omega[d + r][d + r] = value;
				ret.omegafixed[d + r][d + r] = 2;
			}

		/* normal diagonal omega block */
		} else if (type == OMEGA_DIAG) {
			forcount(r, ndim) {
				ret.omega[d + r][d + r] = v[n];
				++n;
			}

		} else
			fatal(0, "invalid OMEGA block type (%i)\n", type);

		/* make sure there are no extra values in list we are ignoring */
		for (int i=n; i<OPENPMX_OMEGABLOCKSIZE_MAX; i++) {
			if (v[i] != 0.)
				fatal(0, "excess value (%f) in omega block\n", v[i]);
		}

		/* set offset for the next block */
		d += ndim; 
	}

	/* find out which elements of omega matrix are fixed */
	forcount(i, ret.nomega) {
		forcount(j, ret.nomega) {
			let v = ret.omega[i][j];
			if (i == j && v <= 0.) {
				ret.omega[i][j] = fabs(v);
				if (ret.omegafixed[i][j] == 2)
					fatal(0, "variances on diagonal of a SAME block cannot be fixed\n");
				ret.omegafixed[i][j] = 1;
			} else if (v == 0.)
				ret.omegafixed[i][j] = 1;
		}
	}

	/* make a copy of the sigma */
	/* we have to check this backwards because there can be zero values at
	 * lower indicies, after that we can go forward */
	for (int i=OPENPMX_SIGMA_MAX-1; i>=0; i--) {
		ret.nsigma = i + 1;
		let v = sigma[i];
		/* test for non-zero, or negative zero */
		if (v != 0.)
			break;
		if (signbit(v))
			break;
	}
	forcount(i, ret.nsigma) {
		let v = sigma[i];

		if (!isfinite(v))
			fatal(0, "invalid sigma %g\n", v);
		
		ret.sigma[i] = fabs(v);
		ret.sigmafixed[i] = (v <= 0.) ? 1 : 0;
	}

	/* initial results are invalid */
	ret.result = (PMXRESULT) { .objfn = DBL_MAX,
							   .type = OBJFN_INVALID,
							   .nparam = 0,
							   .neval = 0 };

	return ret;
}

void extfile_header(FILE * f,
					const POPMODEL* const popmodel,
					const bool _offset1)
{

	char temp[1024];
	let indexoffset = _offset1 ? 1 : 0;

	fprintf(f, OPENPMX_SFORMAT, "ITERATION");
	let ntheta = popmodel->ntheta;
	forcount(i, ntheta) {
		sprintf(temp, "THETA%i", i + indexoffset);
		fprintf(f, OPENPMX_HEADER_FORMAT, temp);
	}
	let nsigma = popmodel->nsigma;
	forcount(i, nsigma) {
		sprintf(temp, "SIGMA(%i,%i)", i + indexoffset, i + indexoffset);
		fprintf(f, OPENPMX_HEADER_FORMAT, temp);
	}
	forcount(i, popmodel->nomega) {
		forcount(j, i+1) {
			sprintf(temp, "OMEGA(%i,%i)", i + indexoffset, j + indexoffset);
			fprintf(f, OPENPMX_HEADER_FORMAT, temp);
		}
	}
	fprintf(f, OPENPMX_HEADER_FORMAT, "OBJ");
	fprintf(f, OPENPMX_SFORMAT, "RUNTIME");
	fprintf(f, OPENPMX_SFORMAT, "INEVAL");
	fprintf(f, "\n");
	fflush(f);
}

void extfile_append(FILE* f, const POPMODEL* const popmodel, const double runtime_s, const int iter, const int ineval)
{
	fprintf(f, OPENPMX_IFORMAT, iter);

	let ntheta = popmodel->ntheta;
	forcount(i, ntheta)
		fprintf(f, OPENPMX_TABLE_FORMAT, popmodel->theta[i]);

	let nsigma = popmodel->nsigma;
	forcount(i, nsigma)
		fprintf(f, OPENPMX_TABLE_FORMAT, popmodel->sigma[i]);

	forcount(i, popmodel->nomega)
		forcount(j, i+1)
			fprintf(f, OPENPMX_TABLE_FORMAT, popmodel->omega[i][j]);

	fprintf(f, OPENPMX_TABLE_FORMAT, popmodel->result.objfn);
	fprintf(f, OPENPMX_TABLE_FORMAT, runtime_s);
	fprintf(f, OPENPMX_IFORMAT, ineval);
	fprintf(f, "\n");
	fflush(f);
}

void extfile_trailer(FILE* f, const POPMODEL* const popmodel, const double runtime_s, const int ineval)
{
	fprintf(f, OPENPMX_IFORMAT, -1000000000);
	let ntheta = popmodel->ntheta;
	forcount(i, ntheta)
		fprintf(f, OPENPMX_TABLE_FORMAT, popmodel->theta[i]);
	let nsigma = popmodel->nsigma;
	forcount(i, nsigma)
		fprintf(f, OPENPMX_TABLE_FORMAT, popmodel->sigma[i]);
	forcount(i, popmodel->nomega)
		forcount(j, i+1)
			fprintf(f, OPENPMX_TABLE_FORMAT, popmodel->omega[i][j]);
	fprintf(f, OPENPMX_TABLE_FORMAT, popmodel->result.objfn);
	fprintf(f, OPENPMX_TABLE_FORMAT, runtime_s);
	fprintf(f, OPENPMX_IFORMAT, ineval);
	fprintf(f, "\n");

	fprintf(f, OPENPMX_IFORMAT, -1000000006);
	forcount(i, ntheta) {
		let est = (popmodel->thetaestim[i] == FIXED) ? 1. : 0.;
		fprintf(f, OPENPMX_TABLE_FORMAT, est);
	}
	forcount(i, nsigma)
		fprintf(f, OPENPMX_TABLE_FORMAT, (double)popmodel->sigmafixed[i]);

	forcount(i, popmodel->nomega)
		forcount(j, i+1)
			fprintf(f, OPENPMX_TABLE_FORMAT, (double)popmodel->omegafixed[i][j]);
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, OPENPMX_IFORMAT, 0);
	fprintf(f, "\n");

	/* NON-STANDARD!!! DIFFERENT THAN NONMEM
	 * add extra line with theta upper and lower limits */
	fprintf(f, OPENPMX_IFORMAT, -2000000001);
	forcount(i, ntheta)
		fprintf(f, OPENPMX_TABLE_FORMAT, popmodel->lower[i]);
	forcount(i, nsigma)
		fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	forcount(i, popmodel->nomega)
		forcount(j, i+1)
			fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, OPENPMX_IFORMAT, 0);
	fprintf(f, "\n");

	fprintf(f, OPENPMX_IFORMAT, -2000000002);
	forcount(i, ntheta)
		fprintf(f, OPENPMX_TABLE_FORMAT, popmodel->upper[i]);
	forcount(i, nsigma)
		fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	forcount(i, popmodel->nomega)
		forcount(j, i+1)
			fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, OPENPMX_IFORMAT, 0);
	fprintf(f, "\n");
	fflush(f);
}

static void info_iteration(FILE* f1,
						   const double runtime_s,
						   const POPMODEL* popmodel,
						   const char* suffix)
{
	let _objfn = popmodel->result.objfn;
	let neval = popmodel->result.neval;
	info(f1, "time %.3f neval %i objfn %.6f", runtime_s, neval, _objfn);
	if (suffix)
		info(f1, "%s", suffix);
	info(f1, "\n");
}

void popmodel_information(FILE* f2, const POPMODEL* const popmodel, const double timestamp)
{
	const char* message = 0;
	switch (popmodel->result.type) {
		case OBJFN_INVALID:		message = "invalid";	break;
		case OBJFN_CURRENT:		message = "current";	break;
		case OBJFN_FINAL:		message = "final";		break;
		case OBJFN_EVALUATE:	message = "evaluate";	break;
		default:
			assert(0);
	}
	assert(message);
	info(f2, "popmodel %s\n", message);

	if (popmodel->result.type == OBJFN_EVALUATE ||
		popmodel->result.type == OBJFN_CURRENT ||
		popmodel->result.type == OBJFN_FINAL) {
		if (timestamp != DBL_MAX)
			info(f2, "time %.3f ", timestamp);
		info(f2, "neval %i ", popmodel->result.neval);
		info(f2, "objfn %.6f ", popmodel->result.objfn);
		info(f2, "nparam %i\n", popmodel->result.nparam);
	}

	/* info about theta */
	let ntheta = popmodel->ntheta;
	info(f2, "THETA(%i)\n", ntheta);
	forcount(i, ntheta) {
		let v = popmodel->theta[i];
		info(f2, OPENPMX_FFORMAT "\n", v);
	}

	/* info about omega matrix */
	let nomega = popmodel->nomega;
	info(f2, "OMEGA(%i)\n", nomega);
	forcount(i, nomega) {
		forcount(j, i+1) 
			info(f2, OPENPMX_FFORMAT " ", popmodel->omega[i][j]);
		info(f2, "\n");
	}

	/* info about sigma */
	let nsigma = popmodel->nsigma;
	info(f2, "SIGMA(%i)\n", nsigma);
	forcount(i, nsigma) 
		info(f2, OPENPMX_FFORMAT " ", popmodel->sigma[i]);
	info(f2, "\n");
}

void popmodel_initcode(FILE* f2, const POPMODEL* const popmodel)
{
	/* info about theta */
	openpmx_printf(f2, 0, 0, "$THETA\n");
	forcount(i, popmodel->ntheta) {
		let l = popmodel->lower[i];
		let v = popmodel->theta[i];
		let u = popmodel->upper[i];
		let s = (popmodel->thetaestim[i] == ESTIMATE) ? "ESTIMATE" : "FIXED";
		openpmx_printf(f2, 0, 0, "\t{" OPENPMX_FFORMAT ",\t" OPENPMX_FFORMAT ",\t" OPENPMX_FFORMAT ",\t %s },\n", l, v, u, s);
	}

	var offset = 0;
	let nblock = popmodel->nblock;
	forcount(k, nblock) {
		let ndim = popmodel->blockdim[k];
		let type = popmodel->blocktype[k];
		switch (type) {

			case OMEGA_DIAG:
				openpmx_printf(f2, 0, 0, "$OMEGA( ");
				forcount(i, ndim) {
					let v = popmodel->omega[offset + i][offset + i];
					let s = popmodel->omegafixed[offset + i][offset + i] ? -1. : 1.;
					openpmx_printf(f2, 0, 0, " %g", v * s);
					if (i != ndim - 1)
						openpmx_printf(f2, 0, 0, ",");
					openpmx_printf(f2, 0, 0, "%s", s < 0. ? "/*FIX*/" :  "");
				}
				openpmx_printf(f2, 0, 0, ")\n");
				break;

			case OMEGA_BLOCK:
				openpmx_printf(f2, 0, 0, "$OMEGABLOCK( ");
				forcount(i, ndim) {
					forcount(j, i+1) {
						let v = popmodel->omega[offset + i][offset + j];
						let s = (i == j && popmodel->omegafixed[offset + i][offset + j]) ? -1. : 1.;
						if (i != 0)
							openpmx_printf(f2, 0, 0, "\t");
						openpmx_printf(f2, 0, 0, "%g", v * s);
						if (i != j || i != ndim - 1)
							openpmx_printf(f2, 0, 0, ",");
						openpmx_printf(f2, 0, 0, "%s", s < 0. ? "/*FIX*/" :  "");
					}
					if (i == ndim - 1)
						openpmx_printf(f2, 0, 0, ")");
				}
				openpmx_printf(f2, 0, 0, "\n");
				break;

			case OMEGA_SAME:
				openpmx_printf(f2, 0, 0, "$OMEGASAME(%i)", ndim);
				break;

			default:
				fatal(f2, "Invalid block type (%i) if size %i\n", type, ndim);
				break;
		}
		offset += ndim;

		if (k != nblock - 1)
			openpmx_printf(f2, 0, 0, "\n");
	}

	/* info about sigma */
	openpmx_printf(f2, 0, 0, "$SIGMA(");
	let nsigma = popmodel->nsigma;
	forcount(i, nsigma) {
		let v = popmodel->sigma[i];
		let s = (popmodel->sigmafixed[i]) ? -1. : 1.;
		openpmx_printf(f2, 0, 0, "\t%g", s * v);
		if (i != nsigma - 1)
			openpmx_printf(f2, 0, 0, ",");
		openpmx_printf(f2, 0, 0, "%s", s < 0. ? "/*FIX*/" :  "");
	}
	openpmx_printf(f2, 0, 0, ")\n");
}

void popmodel_eval_information(const POPMODEL* const popmodel,
							   const double runtime_s,
							   const int neval,
							   const bool details,
							   FILE* outstream,
							   FILE* extstream,
							   const char* suffix)
{
	if (details)
		popmodel_information(outstream, popmodel, runtime_s);

	info_iteration(outstream, runtime_s, popmodel, suffix);

	if (extstream) 
		extfile_append(extstream, popmodel, runtime_s, popmodel->result.neval, neval);
}


