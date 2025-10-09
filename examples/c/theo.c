/*
 * Theophylline oral from Bae and Yim from TCP 2016
 * compile and run with:
 * gcc -W -Wall -Wextra -O2 theo.c -I../../include -I../../src -lgsl -lgslcblas -lm; ./a.out
*/

#include <stdlib.h>
#include <math.h>
#include <stddef.h>

#include "openpmx.h"

typedef struct IMODEL {
	double V, K, KA;
} IMODEL;

typedef struct RECORD {
    double ID, AMT, TIME, DV, BWT;
} RECORD;

static void imodel_init(IMODEL* const _imodel,
						ADVANSTATE* const _advanstate,
						const POPPARAM* const _popparam)
{
	const double* theta = _popparam->theta;
	const double* eta = _popparam->eta;

	const RECORD* record = _advanstate->record;
	const double BWT = record->BWT;

	const double V  = theta[0]*exp(eta[0]) * pow(BWT/70., 1);
    const double K  = theta[1]*exp(eta[1]) * pow(BWT/70., -0.25);
    const double KA = theta[2]*exp(eta[2]) * pow(BWT/70., -0.25);
    
	_imodel-> V = V;
	_imodel-> K = K;
	_imodel-> KA = KA;
}

typedef struct PREDICTVARS {
	/* nothing */
} PREDICTVARS;

static double imodel_predict(const IMODEL* const _imodel,
							 const PREDICTSTATE* const _predictstate,
							 const POPPARAM* const _popparam,
							 const double* const err,
							 PREDICTVARS* _predparams)
{
	(void)_popparam;
	(void)_predparams;
	
	const RECORD* record = _predictstate->record;
	const double TIME = record->TIME;

	const double V = _imodel->V;
	const double K = _imodel->K;
	const double KA = _imodel->KA;

	const double DOSE=320.;
    const double IPRED = DOSE / V * KA / (KA - K) * (exp(-K * TIME)-exp(-KA * TIME));
    
    return IPRED * (1 + err[0]) + err[1];
}

#define ARRAYSIZE(a) (sizeof(a)/sizeof(a[0]))

extern RECORD data[132]; /* forward declaration */

int main(void)
{
	OPENPMX openpmx = (OPENPMX) {
		.filename = "theo.c",
		.data = {
			.writeable = 0,
			.records = data,
			.nrecords = ARRAYSIZE(data),
			.recordfields = {
				.size = sizeof(RECORD),
				.field = { 
					{ .name="ID", .offset=offsetof(RECORD, ID) },
					{ .name="AMT", .offset=offsetof(RECORD, AMT) },
					{ .name="TIME", .offset=offsetof(RECORD, TIME) },
					{ .name="DV", .offset=offsetof(RECORD, DV) },
					{ .name="BWT", .offset=offsetof(RECORD, BWT) },
				},
			},
		},
		.advan = {
			.init = imodel_init,
			.predict = imodel_predict,
			.diffeqn = 0,
			.imodelfields = {
				.size = sizeof(IMODEL),
				.field = {
					{ .name="V", .offset = offsetof(IMODEL, V) },
					{ .name="K", .offset = offsetof(IMODEL, K) },
					{ .name="KA", .offset = offsetof(IMODEL, KA) },
				},
			},
			.predictfields = {
				.size = sizeof(PREDICTVARS),
				.field = { },
			},
			.method = pmx_advan_pred,
			.firstonly = true,
		},
		.theta = {	{ 0.1, 15, 50, ESTIMATE },
					{ 0.0001, 0.1, 1, ESTIMATE },
					{ 0.0001, 0.5, 5, ESTIMATE },
		},
		.omega = {	{ OMEGA_BLOCK, 3, { .1,
										0.01, .1,
										0.01, 0.01, .1 } },
		},
		.sigma = { 1, 1 },
	};

	pmx_estimate(&openpmx, &(ESTIMCONFIG) { });
	pmx_cleanup(&openpmx);

	return EXIT_SUCCESS;
}

/* include sources directly instead of linking to library */
#include "dataconfig/dataconfig.c"
#include "advan/advan.c"
#include "advan/pred.c"
#include "bobyqa/bobyqa.c"
#include "pmxstate.c"
#include "idata.c"
#include "print.c"
#include "options.c"
#include "popmodel.c"
#include "omegainfo.c"
#include "linalg.c"
#include "ievaluate.c"
#include "stage1.c"
#include "scatter.c"
#include "estimate.c"
#include "checkout.c"
#include "encode.c"
#include "utils/vector.c"
#include "utils/various.c"

RECORD data[132] = {
	{	1,	4.02,	0.,	.74,	79.6,	},
	{	1,	0.,	0.25,	2.84,	0.,	},
	{	1,	0.,	0.57,	6.57,	0.,	},
	{	1,	0.,	1.12,	10.5,	0.,	},
	{	1,	0.,	2.02,	9.66,	0.,	},
	{	1,	0.,	3.82,	8.58,	0.,	},
	{	1,	0.,	5.1,	8.36,	0.,	},
	{	1,	0.,	7.03,	7.47,	0.,	},
	{	1,	0.,	9.05,	6.89,	0.,	},
	{	1,	0.,	12.12,	5.94,	0.,	},
	{	1,	0.,	24.37,	3.28,	0.,	},
	{	2,	4.4,	0.,	0.,	72.4,	},
	{	2,	0.,	.27,	1.72,	0.,	},
	{	2,	0.,	.52,	7.91,	0.,	},
	{	2,	0.,	1.,	8.31,	0.,	},
	{	2,	0.,	1.92,	8.33,	0.,	},
	{	2,	0.,	3.5,	6.85,	0.,	},
	{	2,	0.,	5.02,	6.08,	0.,	},
	{	2,	0.,	7.03,	5.4,	0.,	},
	{	2,	0.,	9.,	4.55,	0.,	},
	{	2,	0.,	12.,	3.01,	0.,	},
	{	2,	0.,	24.3,	.90,	0.,	},
	{	3,	4.53,	0.,	0.,	70.5,	},
	{	3,	0.,	.27,	4.4,	0.,	},
	{	3,	0.,	.58,	6.9,	0.,	},
	{	3,	0.,	1.02,	8.2,	0.,	},
	{	3,	0.,	2.02,	7.8,	0.,	},
	{	3,	0.,	3.62,	7.5,	0.,	},
	{	3,	0.,	5.08,	6.2,	0.,	},
	{	3,	0.,	7.07,	5.3,	0.,	},
	{	3,	0.,	9.,	4.9,	0.,	},
	{	3,	0.,	12.15,	3.7,	0.,	},
	{	3,	0.,	24.17,	1.05,	0.,	},
	{	4,	4.4,	0.,	0.,	72.7,	},
	{	4,	0.,	.35,	1.89,	0.,	},
	{	4,	0.,	.6,	4.6,	0.,	},
	{	4,	0.,	1.07,	8.6,	0.,	},
	{	4,	0.,	2.13,	8.38,	0.,	},
	{	4,	0.,	3.5,	7.54,	0.,	},
	{	4,	0.,	5.02,	6.88,	0.,	},
	{	4,	0.,	7.02,	5.78,	0.,	},
	{	4,	0.,	9.02,	5.33,	0.,	},
	{	4,	0.,	11.98,	4.19,	0.,	},
	{	4,	0.,	24.65,	1.15,	0.,	},
	{	5,	5.86,	0.,	0.,	54.6,	},
	{	5,	0.,	.3,	2.02,	0.,	},
	{	5,	0.,	.52,	5.63,	0.,	},
	{	5,	0.,	1.,	11.4,	0.,	},
	{	5,	0.,	2.02,	9.33,	0.,	},
	{	5,	0.,	3.5,	8.74,	0.,	},
	{	5,	0.,	5.02,	7.56,	0.,	},
	{	5,	0.,	7.02,	7.09,	0.,	},
	{	5,	0.,	9.1,	5.9,	0.,	},
	{	5,	0.,	12.,	4.37,	0.,	},
	{	5,	0.,	24.35,	1.57,	0.,	},
	{	6,	4.,	0.,	0.,	80.,	},
	{	6,	0.,	.27,	1.29,	0.,	},
	{	6,	0.,	.58,	3.08,	0.,	},
	{	6,	0.,	1.15,	6.44,	0.,	},
	{	6,	0.,	2.03,	6.32,	0.,	},
	{	6,	0.,	3.57,	5.53,	0.,	},
	{	6,	0.,	5.,	4.94,	0.,	},
	{	6,	0.,	7.,	4.02,	0.,	},
	{	6,	0.,	9.22,	3.46,	0.,	},
	{	6,	0.,	12.1,	2.78,	0.,	},
	{	6,	0.,	23.85,	.92,	0.,	},
	{	7,	4.95,	0.,	.15,	64.6,	},
	{	7,	0.,	.25,	.85,	0.,	},
	{	7,	0.,	.5,	2.35,	0.,	},
	{	7,	0.,	1.02,	5.02,	0.,	},
	{	7,	0.,	2.02,	6.58,	0.,	},
	{	7,	0.,	3.48,	7.09,	0.,	},
	{	7,	0.,	5.,	6.66,	0.,	},
	{	7,	0.,	6.98,	5.25,	0.,	},
	{	7,	0.,	9.,	4.39,	0.,	},
	{	7,	0.,	12.05,	3.53,	0.,	},
	{	7,	0.,	24.22,	1.15,	0.,	},
	{	8,	4.53,	0.,	0.,	70.5,	},
	{	8,	0.,	.25,	3.05,	0.,	},
	{	8,	0.,	0.52,	3.05,	0.,	},
	{	8,	0.,	.98,	7.31,	0.,	},
	{	8,	0.,	2.02,	7.56,	0.,	},
	{	8,	0.,	3.53,	6.59,	0.,	},
	{	8,	0.,	5.05,	5.88,	0.,	},
	{	8,	0.,	7.15,	4.73,	0.,	},
	{	8,	0.,	9.07,	4.57,	0.,	},
	{	8,	0.,	12.1,	3.,	0.,	},
	{	8,	0.,	24.12,	1.25,	0.,	},
	{	9,	3.1,	.0,	.0,	86.4,	},
	{	9,	0.,	.3,	7.37,	0.,	},
	{	9,	0.,	.63,	9.03,	0.,	},
	{	9,	0.,	1.05,	7.14,	0.,	},
	{	9,	0.,	2.02,	6.33,	0.,	},
	{	9,	0.,	3.53,	5.66,	0.,	},
	{	9,	0.,	5.02,	5.67,	0.,	},
	{	9,	0.,	7.17,	4.24,	0.,	},
	{	9,	0.,	8.8,	4.11,	0.,	},
	{	9,	0.,	11.6,	3.16,	0.,	},
	{	9,	0.,	24.43,	1.12,	0.,	},
	{	10,	5.5,	0.,	.24,	58.2,	},
	{	10,	0.,	.37,	2.89,	0.,	},
	{	10,	0.,	.77,	5.22,	0.,	},
	{	10,	0.,	1.02,	6.41,	0.,	},
	{	10,	0.,	2.05,	7.83,	0.,	},
	{	10,	0.,	3.55,	10.21,	0.,	},
	{	10,	0.,	5.05,	9.18,	0.,	},
	{	10,	0.,	7.08,	8.02,	0.,	},
	{	10,	0.,	9.38,	7.14,	0.,	},
	{	10,	0.,	12.1,	5.68,	0.,	},
	{	10,	0.,	23.7,	2.42,	0.,	},
	{	11,	4.92,	0.,	0.,	65.,	},
	{	11,	0.,	.25,	4.86,	0.,	},
	{	11,	0.,	.5,	7.24,	0.,	},
	{	11,	0.,	.98,	8.,	0.,	},
	{	11,	0.,	1.98,	6.81,	0.,	},
	{	11,	0.,	3.6,	5.87,	0.,	},
	{	11,	0.,	5.02,	5.22,	0.,	},
	{	11,	0.,	7.03,	4.45,	0.,	},
	{	11,	0.,	9.03,	3.62,	0.,	},
	{	11,	0.,	12.12,	2.69,	0.,	},
	{	11,	0.,	24.08,	.86,	0.,	},
	{	12,	5.3,	0.,	0.,	60.5,	},
	{	12,	0.,	.25,	1.25,	0.,	},
	{	12,	0.,	.5,	3.96,	0.,	},
	{	12,	0.,	1.,	7.82,	0.,	},
	{	12,	0.,	2.,	9.72,	0.,	},
	{	12,	0.,	3.52,	9.75,	0.,	},
	{	12,	0.,	5.07,	8.57,	0.,	},
	{	12,	0.,	7.07,	6.59,	0.,	},
	{	12,	0.,	9.03,	6.11,	0.,	},
	{	12,	0.,	12.05,	4.57,	0.,	},
	{	12,	0.,	24.15,	1.17,	0.,	},
};

