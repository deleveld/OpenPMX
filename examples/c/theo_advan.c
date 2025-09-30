/*
 * Theophylline oral from Bae and Yim from TCP 2016
 * compile and run with:
 * gcc -W -Wall -Wextra -O2 theo_advan.c -I../../include -I../../src -lgsl -lgslcblas -lm; ./a.out
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
	double IPRED;
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
    _predparams->IPRED = IPRED;
    
    return IPRED * (1 + err[0]) + err[1];
}

#define ARRAYSIZE(a) (sizeof(a)/sizeof(a[0]))

extern RECORD data[11]; /* forward declaration */

/* include sources directly instead of linking to library */
#include "dataconfig/dataconfig.c"
#include "dataconfig/recordinfo.c"
#include "advan/advan.c"
#include "advan/pred.c"
#include "print.c"
#include "popmodel.c"
#include "ievaluate.h"
#include "utils/vector.c"

int main(void)
{
	OPENPMX openpmx = (OPENPMX) {
		.data = {
			.writeable = data,
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
			.imodelfields = { },
			.predictfields = { },
			.method = pmx_advan_pred,
			.firstonly = true,
		},
		.theta = {	{          0.1,	      32.7707,	           50,	 ESTIMATE },
					{       0.0001,	    0.0867683,	            1,	 ESTIMATE },
					{       0.0001,	      1.51591,	            5,	 ESTIMATE },
		},
		.omega = {	{ OMEGA_BLOCK, 3, { 0.0108327,
										0.0142938,	0.018945,
										0.00866303,	0.00981528,	0.454379 } },
		},
		.sigma = { 0.0168844,	0.0732164 },
	};

	/* get advancer function table, memory needed to iterate, and construct advancer */
	let advanfuncs = advanfuncs_alloc(&openpmx.data, &openpmx.advan);
	ADVAN* advan = calloc(advanfuncs->advan_size, 1);	
	advanfuncs->construct(advan, advanfuncs);
	
	/* arguments needed to iterate */
	double eta[OPENPMX_OMEGA_MAX] = { };
	let popmodel = popmodel_init(&openpmx);
	let ievaluate_args = (IEVALUATE_ARGS) {
		.record = data,
		.nrecord = ARRAYSIZE(data),
		.advanfuncs = advanfuncs,
		.popparam = { 
			.theta = popmodel.theta,
			.ntheta = popmodel.ntheta,
			.eta = eta,
			.nomega = popmodel.nomega,
			.sigma = popmodel.sigma,
			.nsigma = popmodel.nsigma,
			.nstate = advanfuncs->nstate,
		},
		.logstream = 0,
	};
	
	/* scratch memory needed to iterate */
	double errarray[OPENPMX_SIGMA_MAX] = { };
	IMODEL imodel;
	PREDICTVARS predictvars;

	/* advance over individuals data */
	let predict = openpmx.advan.predict;
	forcount(i, ARRAYSIZE(data)) {
		let ptr = &data[i];
		let predictstate = advan_advance(advan, &imodel, ptr, &ievaluate_args.popparam);
		let yhat = predict(&imodel, &predictstate, &ievaluate_args.popparam, errarray, &predictvars);
		
		printf("%f %f %f %f %f\n", ptr->ID, ptr->TIME, ptr->DV, yhat, predictvars.IPRED);
	}
	advanfuncs->destruct(advan);
	free(advan);

	return EXIT_SUCCESS;
}

RECORD data[11] = {
	{	1,	4.02,	0.,	.74,	179.6,	},
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
};

