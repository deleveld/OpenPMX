/* 
 * This file is part of OpenPMX (https://github.com/deleveld/openpmx).
 * Copyright (c) 2022 Douglas Eleveld.
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
#include <values.h>

#include "popmodel.h"
#include "print.h"
#include "defines.h"
#include "utils/c22.h"
#include "utils/various.h"

#include "openpmx_compile_options.h"

/* TODO: why does this give a com piler error */
POPMODEL popmodel_init(const THETA theta[static OPENPMX_THETA_MAX],
					   const OMEGA* omegablocks, // not sure why static here causes warning [static OPENPMX_OMEGABLOCK_MAX],
					   const double sigma[static OPENPMX_SIGMA_MAX])
{
	POPMODEL ret = { 0 };
	ret.ntheta = 0;
	ret.nblock = 0;
	ret.nomega = 0;
	ret.nsigma = 0;

	/* initialize to invalid values to catch errors */
	forcount(i, OPENPMX_THETA_MAX) {
		ret.lower[i] = DBL_MAX;
		ret.theta[i] = DBL_MAX;
		ret.upper[i] = DBL_MAX;
		ret.thetaestim[i] = FIXED;
	}
	forcount(i, OPENPMX_OMEGABLOCK_MAX) {
		ret.blocktype[i] = OMEGA_INVALID;
		ret.blockdim[i] = 0;
	}
	forcount(i, OPENPMX_SIGMA_MAX) {
		ret.sigma[i] = DBL_MAX;
		ret.sigmafixed[i] = 1;
	}

	/* count thetas, omegas, sigmas */
	forcount(i, OPENPMX_THETA_MAX) {
		let est = theta[i].type;
		if (est == THETA_INVALID) 
			break;
		ret.ntheta = i + 1;
		
		ret.lower[i] = theta[i].lower;
		ret.theta[i] = theta[i].value;
		ret.upper[i] = theta[i].upper;

		if (theta[i].value < theta[i].lower ||
			theta[i].value > theta[i].upper)
			fatal(0, "THETA outside boundary { %g, %g, %g }\n", theta[i].lower, theta[i].value, theta[i].upper);

		if (est == ESTIMATE) {
			if (theta[i].value == theta[i].lower ||
				theta[i].value == theta[i].upper)
				fatal(0, "THETA %i at boundary { %g, %g, %g, ESTIMATE }\n", i, theta[i].lower, theta[i].value, theta[i].upper);
		}
		ret.thetaestim[i] = est;
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
		fatal(0, "Omega size (%i) is too large, Max is %i\n", ret.nomega, OPENPMX_OMEGA_MAX);

	/* initialize the omega matrix from the blocks */
	var d = 0;
	forcount(i, ret.nblock) {
		let ndim = omegablocks[i].ndim;
		let v = omegablocks[i].values;
		let type = omegablocks[i].type;
		if (type == OMEGA_BLOCK) {
			var n = 0;
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
		} else if (type == OMEGA_SAME) {
			/* get entries from elsewhere set omega elements */
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
		} else if (type == OMEGA_DIAG) {
			/* normal diagonal omega block */
			var n = 0;
			forcount(r, ndim) {
				ret.omega[d + r][d + r] = v[n];
				++n;
			}
		} else
			fatal(0, "illegal OMEGATYPE (%i)\n", type);

		d += ndim;
	}

	/* find out which elements of omega matrix are fixed */
	forcount(i, ret.nomega) {
		forcount(j, ret.nomega) {
			let v = ret.omega[i][j];
			if (i == j && v <= 0.) {
				ret.omega[i][j] = fabs(v);
				if (ret.omegafixed[i][j] == 2)
					fatal(0, "variances on diagonal that are set as fixed cannot be part of a SAME block\n");
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
		let fabsv = fabs(v);
		if (fabsv != 0.)
			break;
	}
	forcount(i, ret.nsigma) {
		let v = sigma[i];
		let fabsv = fabs(v);
		ret.sigma[i] = fabsv;
		ret.sigmafixed[i] = (v <= 0.) ? 1 : 0;
	}

	/* initial results are invalid */
	ret.result = (PMXRESULT) { .objfn = DBL_MAX,
							   .type = OBJFN_INVALID,
							   .nfunc = 0 };

	return ret;
}

void extfile_header(const char* filename,
					const POPMODEL* const popmodel,
					const bool offset_1)
{
	var f = results_fopen(filename, OPENPMX_EXTFILE, "w");
	assert(f);

	char temp[1024];
	let indexoffset = offset_1 ? 1 : 0;

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
	fprintf(f, OPENPMX_HEADER_FORMAT, "MAXD");
	fprintf(f, "\n");

	fclose(f);
}

void extfile_append(const char* filename,
					const POPMODEL* const popmodel,
					const double maxd)
{
	var f = results_fopen(filename, OPENPMX_EXTFILE, "a");
	if (f == 0)
		fatal(0, "%s: could not open file %s extension %s\n", __func__, filename, OPENPMX_EXTFILE);

	let iter = popmodel->result.nfunc;
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
	fprintf(f, OPENPMX_TABLE_FORMAT, maxd);
	fprintf(f, "\n");

	fclose(f);
}

void extfile_trailer(const char* filename, const POPMODEL* const popmodel)
{
	var f = results_fopen(filename, OPENPMX_EXTFILE, "a");
	assert(f);

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

	/* objfn and maxd */
	fprintf(f, "  ");
	fprintf(f, OPENPMX_TABLE_FORMAT, popmodel->result.objfn);
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
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

	/* objfn and maxd */
	fprintf(f, "  ");
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
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
	fprintf(f, "  ");
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, "\n");

	fprintf(f, OPENPMX_IFORMAT, -2000000002);
	forcount(i, ntheta)
		fprintf(f, OPENPMX_TABLE_FORMAT, popmodel->upper[i]);
	forcount(i, nsigma)
		fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	forcount(i, popmodel->nomega)
		forcount(j, i+1)
			fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, "  ");
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, OPENPMX_TABLE_FORMAT, 0.);
	fprintf(f, "\n");

	fclose(f);
}

void info_iteration(FILE* f1,
					const double runtime_s,
					const double d1,
					const POPMODEL* popmodel)
{
	let _objfn = popmodel->result.objfn;
	let nfunc = popmodel->result.nfunc;
	info(f1, "time: %.3f nfunc: %i objfn: %.6f", runtime_s, nfunc, _objfn);
	if (d1 != 0.)
		info(f1, " d: %g", d1);
	info(f1, "\n");
}

void print_iteration(FILE* f1,
					 FILE* f2,
					 const POPMODEL* popmodel,
					 const int xlength,
					 const double* const x)
{
	openpmx_printf(f1, f2, 0, "param:");

	let etheta = popmodel->theta;
	let ntheta = popmodel->ntheta;
	let thetaestim = popmodel->thetaestim;
	forcount(j, ntheta) {
		if (thetaestim[j] != FIXED) {
			let v = etheta[j];
			openpmx_printf(f1, f2, 0, " % -11.4e", v);
		}
	}
	let sigma = popmodel->sigma;
	let nsigma = popmodel->nsigma;
	let sigmafixed = popmodel->sigmafixed;
	forcount(j, nsigma) {
		if (sigmafixed[j] == 0) {
			let v = sigma[j];
			openpmx_printf(f1, f2, 0, " % -11.4e", v);
		}
	}
	let nomega = popmodel->nomega;
	forcount(j, nomega) {
		for (var k=0; k<=j; k++) {
			let v = popmodel->omega[j][k];
			let f = popmodel->omegafixed[j][k];
			if (f == 0 && v != 0)
				openpmx_printf(f1, f2, 0, " % -11.4e", v);
		}
	}
	if (x && xlength > 0) {
		openpmx_printf(f1, f2, 0, "\ntpara:");
		forcount(j, xlength)
			openpmx_printf(f1, f2, 0, " % -11.4e", x[j]);
		openpmx_printf(f1, f2, 0, "\nd    :");
		forcount(j, xlength)
			openpmx_printf(f1, f2, 0, " % -11.4e", fabs(x[j]));
	}
	openpmx_printf(f1, f2, 0, "\n");
}

void iterfile_popmodel_information(FILE* f2, const POPMODEL* const popmodel)
{
	if (popmodel->result.type == OBJFN_EVALUATE ||
		popmodel->result.type == OBJFN_CURRENT ||
		popmodel->result.type == OBJFN_FINAL) {
		info(f2, "nfunc: %i\n", popmodel->result.nfunc);
		info(f2, "objfn: %.6f\n", popmodel->result.objfn);
	}

	const char* message = 0;
	switch (popmodel->result.type) {
		case OBJFN_INVALID:		message = "invalid:\n";	break;
		case OBJFN_CURRENT:		message = "current:\n";	break;
		case OBJFN_FINAL:		message = "final:\n";	break;
		case OBJFN_PREDICT:		message = "predict:\n";	break;
		case OBJFN_EVALUATE:	message = "evaluate:\n";	break;
		case OBJFN_INITIAL:		message = "initial:\n";		break;
	}
	if (message)
		info(f2, message);

	/* info about theta */
	let ntheta = popmodel->ntheta;
	info(f2, "$THETA\n");
	forcount(i, ntheta) {
		let l = popmodel->lower[i];
		let v = popmodel->theta[i];
		let u = popmodel->upper[i];
		let s = (popmodel->thetaestim[i] == ESTIMATE) ? "ESTIMATE" : "FIXED";
		info(f2, "\t{" OPENPMX_FFORMAT ",\t" OPENPMX_FFORMAT ",\t" OPENPMX_FFORMAT ",\t %s },\n", l, v, u, s);
	}

	/* TODO: This has to be done correctly, in the same blocks the user made
	 * otherwise SAME blocks are screwed up. FIXME!!! */
	/* info about omega matrix */
	let nomega = popmodel->nomega;
	info(f2, "$OMEGABLOCK( ");
	forcount(i, nomega) {
		forcount(j, i+1) {
			let v = popmodel->omega[i][j];
			let s = (i == j && popmodel->omegafixed[i][j]) ? -1. : 1.;
			if (i != 0)
				info(f2, "\t");
			info(f2, "%g", v * s);
			if (i != j || i != nomega - 1)
				info(f2, ",");
			info(f2, "%s", s < 0. ? "/*FIX*/" :  "");
		}
		if (i == nomega - 1)
			info(f2, ")");
		info(f2, "\n");
	}

	/* info about sigma */
	let nsigma = popmodel->nsigma;
	info(f2, "$SIGMA(");
	forcount(i, nsigma) {
		let v = popmodel->sigma[i];
		let s = (popmodel->sigmafixed[i]) ? -1. : 1.;
		info(f2, "\t%g", s * v);
		if (i != nsigma - 1)
			info(f2, ",");
		info(f2, "%s", s < 0. ? "/*FIX*/" :  "");
	}
	info(f2, ")\n");
}

