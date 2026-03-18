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
#include "defines.h"
#include "print.h"
#include "utils/c22.h"
#include "utils/various.h"
#include "utils/errctx.h"

#include <gsl/gsl_math.h>

POPMODEL popmodel_init(const OPENPMX* const pmx, ERRCTX* errctx)
{
	let theta = pmx->theta;
	let omegablocks = pmx->omega;
	let sigma = pmx->sigma;

	POPMODEL ret = { 0 };
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
	forcount(i, OPENPMX_OMEGA_MAX) {
		forcount(j, OPENPMX_OMEGA_MAX) {
			ret.omega[i][j] = NAN;
			ret.omegafixed[i][j] = 0;
		}
	}
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

		if (!gsl_finite(theta[i].lower) ||
			!gsl_finite(theta[i].value) ||
			!gsl_finite(theta[i].upper)) {
			add_errctx(errctx, "THETA invalid { %g, %g, %g }\n", theta[i].lower, theta[i].value, theta[i].upper);
			goto failed;
		}
		
		ret.lower[i] = theta[i].lower;
		ret.theta[i] = theta[i].value;
		ret.upper[i] = theta[i].upper;

		if (est == ESTIMATE) {
			if (theta[i].value < theta[i].lower ||
				theta[i].value > theta[i].upper) {
				add_errctx(errctx, "THETA outside boundary { %g, %g, %g }\n", theta[i].lower, theta[i].value, theta[i].upper);
				goto failed;
			}

			if (theta[i].value == theta[i].lower ||
				theta[i].value == theta[i].upper) {
				add_errctx(errctx, "THETA at boundary { %g, %g, %g, ESTIMATE }\n", theta[i].lower, theta[i].value, theta[i].upper);
				goto failed;
			}
		}
		ret.thetaestim[i] = est;
	}
	for(int i=ret.ntheta; i<OPENPMX_THETA_MAX; i++) {
		if (theta[i].lower != 0. ||
			theta[i].value != 0. ||
			theta[i].upper != 0. ||
			theta[i].type != THETA_INVALID) {
			add_errctx(errctx, "THETA not fully initialized { %g, %g, %g, ??? }\n", theta[i].lower, theta[i].value, theta[i].upper);
			goto failed;
		}
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
	if (ret.nomega > OPENPMX_OMEGA_MAX) {
		add_errctx(errctx, "Omega size (%i) is too large, max is %i\n", ret.nomega, OPENPMX_OMEGA_MAX);
		goto failed;
	}

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

		} else {
			add_errctx(errctx, "invalid OMEGA block type (%i)\n", type);
			goto failed;
		}

		/* make sure there are no extra values in list we are ignoring */
		for (int i=n; i<OPENPMX_OMEGABLOCKSIZE_MAX; i++) {
			if (v[i] != 0.) {
				add_errctx(errctx, "excess value (%f) in omega block\n", v[i]);
				goto failed;
			}
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
				if (ret.omegafixed[i][j] == 2) {
					add_errctx(errctx, "variances on diagonal of a SAME block cannot be fixed\n");
					goto failed;
				}
				ret.omegafixed[i][j] = 1;
			} else if (v == 0.)
				ret.omegafixed[i][j] = 1;
		}
	}

	/* make a copy of the sigma */
	/* last non-zero sigma is the last one */
	ret.nsigma = 0;
	forcount(i, OPENPMX_SIGMA_MAX) {
		let v = sigma[i];
		if (!isnan(v) && v != 0.)
			ret.nsigma = i + 1;
	}
	forcount(i, ret.nsigma) {
		let v = sigma[i];

		if (!gsl_finite(v)) {
			add_errctx(errctx, "invalid sigma[%i] %g\n", i, v);
			goto failed;
		}
		
		ret.sigma[i] = fabs(v);
		ret.sigmafixed[i] = (v <= 0.) ? 1 : 0;
	}

	/* success path */
	/* initial results are invalid */
	ret.result = (PMXRESULT) { 0 };
	return ret;

	/* failure path */
failed:
	ret = (POPMODEL) { 0 };
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
	fprintf(f, OPENPMX_HEADER_FORMAT, "RUNTIME");
	fprintf(f, OPENPMX_SFORMAT, "INEVAL");
	fprintf(f, "\n");

	/* NON-STANDARD!!! DIFFERENT THAN NONMEM
	 * add extra line with theta upper and lower limits */

	fprintf(f, OPENPMX_IFORMAT, OPENPMX_EXTFILE_LOWER_BOUNDS);
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

	fprintf(f, OPENPMX_IFORMAT, OPENPMX_EXTFILE_UPPER_BOUNDS);
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

	fprintf(f, OPENPMX_IFORMAT, OPENPMX_EXTFILE_EXTRA_FIXED);
	forcount(i, ntheta) {
		let est = (popmodel->thetaestim[i] == FIXED) ? 1. : 0.;
		fprintf(f, OPENPMX_TABLE_FORMAT, est);
	}
	forcount(i, nsigma)
		fprintf(f, OPENPMX_TABLE_FORMAT, (double)popmodel->sigmafixed[i]);

	forcount(i, popmodel->nomega) {
		forcount(j, i+1) {
			fprintf(f, OPENPMX_TABLE_FORMAT, (double)popmodel->omegafixed[i][j]);
		}
	}
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, OPENPMX_IFORMAT, 0);
	fprintf(f, "\n");

	fprintf(f, OPENPMX_IFORMAT, OPENPMX_EXTFILE_OMEGA_BLOCKS);
	forcount(i, ntheta) 
		fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	forcount(i, nsigma)
		fprintf(f, OPENPMX_TABLE_FORMAT, 0.);

	int omegablock[OPENPMX_OMEGA_MAX][OPENPMX_OMEGA_MAX] = { 0 };
	var d = 0;
	forcount(i, popmodel->nblock) {
		let ndim = popmodel->blockdim[i];
		let type = popmodel->blocktype[i];
		var n = 0;
		if (type == OMEGA_BLOCK) {
			forcount(r, ndim) {
				/* off diagonal */
				for (var c=0; c<r; c++) {
					omegablock[d + c][d + r] = i + 1;
					omegablock[d + r][d + c] = i + 1;
					++n;
				}
				/* on diagonal */
				omegablock[d + r][d + r] = i + 1;
				++n;
			}

		/* get entries from elsewhere set omega elements */
		} else if (type == OMEGA_SAME) {
			forcount(r, ndim) {
				/* off diagonal */
				for (var c=0; c<r; c++) {
					omegablock[d + c][d + r] = i + 1;
					omegablock[d + r][d + c] = i + 1;
					++n;
				}
				/* on diagonal */
				omegablock[d + r][d + r] = i + 1;
				++n;
			}

		/* normal diagonal omega block */
		} else if (type == OMEGA_DIAG) {
			forcount(r, ndim) {
				omegablock[d + r][d + r] = i + 1;
				++n;
			}
			
		} else
			fatal(0, "invalid OMEGA block type (%i)\n", type);

		/* set offset for the next block */
		d += ndim; 
	}
	forcount(i, popmodel->nomega)
		forcount(j, i+1)
			fprintf(f, OPENPMX_TABLE_FORMAT, (double)omegablock[i][j]);
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, OPENPMX_IFORMAT, 0);
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
#define NONMEM_FINAL_ESTIMATE	-1000000000
#define NONMEM_FIXED_FLAGS		-1000000006

	fprintf(f, OPENPMX_IFORMAT, NONMEM_FINAL_ESTIMATE);
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

	fprintf(f, OPENPMX_IFORMAT, NONMEM_FIXED_FLAGS);
	forcount(i, ntheta) {
		let est = (popmodel->thetaestim[i] == FIXED) ? 1. : 0.;
		fprintf(f, OPENPMX_TABLE_FORMAT, est);
	}
	forcount(i, nsigma)
		fprintf(f, OPENPMX_TABLE_FORMAT, (double)popmodel->sigmafixed[i]);

	forcount(i, popmodel->nomega)
		forcount(j, i+1)
			fprintf(f, OPENPMX_TABLE_FORMAT, popmodel->omegafixed[i][j] ? 1. : 0.);
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
	if (popmodel->result.type != OBJFN_INVALID) {
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
	}
	
	/* info about theta */
	let ntheta = popmodel->ntheta;
	info(f2, "THETA(%i)\n", ntheta);
	forcount(i, ntheta) {
		let v = popmodel->theta[i];
		let f = popmodel->thetaestim[i];
		info(f2, OPENPMX_FFORMAT "%s\n", 
			v, 
			f == FIXED ? "*" : "");
	}

	/* info about omega matrix */
	let nomega = popmodel->nomega;
	info(f2, "OMEGA(%i)\n", nomega);
	forcount(i, nomega) {
		forcount(j, i+1) {
			let v = popmodel->omega[i][j];
			let f = popmodel->omegafixed[i][j];
			info(f2, OPENPMX_FFORMAT "%s", v, f ? "*" : "");
		}
		info(f2, "\n");
	}

	/* info about sigma */
	let nsigma = popmodel->nsigma;
	info(f2, "SIGMA(%i)\n", nsigma);
	forcount(i, nsigma) {
		let v = popmodel->sigma[i];
		let f = popmodel->sigmafixed[i];
		info(f2, OPENPMX_FFORMAT " %s\n", v, f ? "*" : "");
	}
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


