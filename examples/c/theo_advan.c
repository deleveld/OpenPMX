/*
 * Theophylline oral from Bae and Yim from TCP 2016
 * Compile and run with:
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
						ADVANSTATE* const _advanstate)
{
	const POPPARAM* _popparam = _advanstate->current.popparam;
	const double* theta = _popparam->theta;
	const double* eta = _popparam->eta;
	const RECORD* record = _advanstate->current.record;
	const double BWT = record->BWT;

	_imodel->V  = theta[0]*exp(eta[0]) * pow(BWT/70., 1);
    _imodel->K  = theta[1]*exp(eta[1]) * pow(BWT/70., -0.25);
    _imodel->KA = theta[2]*exp(eta[2]) * pow(BWT/70., -0.25);
}

typedef struct PREDICTVARS {
	double IPRED;
} PREDICTVARS;

static double imodel_predict(const IMODEL* const _imodel,
							 const PREDICTSTATE* const _predictstate,
							 const double* const err,
							 PREDICTVARS* _predparams)
{
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

const RECORD data[] = {
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

/* include sources directly instead of linking to library */
#include "dataconfig/dataconfig.c"
#include "advan/advan.c"
#include "advan/pred.c"
#include "utils/vector.c"

int main(void)
{
	OPENPMX openpmx = (OPENPMX) {
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
			.method = pmx_advan_pred,
			.firstonly = true,
		},
		.theta = {	{          0.1,	      32.7707,	           50,	 ESTIMATE },
					{       0.0001,	    0.0867683,	            1,	 ESTIMATE },
					{       0.0001,	      1.51591,	            5,	 ESTIMATE },
		},
		.omega = {	
			{ OMEGA_DIAG, 3, { 0, 0, 0 } }, 
		},
		.sigma = { 0, 0 },
	};

	/* advanfuncs and arguments needed to iterate */
	let advanfuncs = advanfuncs_alloc(&openpmx.data, &openpmx.advan);
	const double theta[] = { 
		openpmx.theta[0].value,
		openpmx.theta[1].value,
		openpmx.theta[2].value,
	};
	double eta[OPENPMX_OMEGA_MAX] = { };
	let popparam = (POPPARAM) { 
		.theta = theta,
		.ntheta = ARRAYSIZE(theta),
		.eta = eta,
		.nomega = 0,
		.sigma = 0,
		.nsigma = 0,
		.nstate = advanfuncs->nstate,
	};
	
	/* get memory needed to iterate, and construct advancer */
	ADVAN* advan = calloc(advanfuncs->advan_size, 1);	
	advanfuncs->construct(advan, advanfuncs);
	/* scratch memory needed to iterate */
	double errarray[OPENPMX_SIGMA_MAX] = { };
	IMODEL imodel;
	PREDICTVARS predictvars;

	/* advance over individuals data */
	let predict = openpmx.advan.predict;
	forcount(i, ARRAYSIZE(data)) {
		let ptr = &data[i];
		let predictstate = advan_advance(advan, &imodel, ptr, &popparam);
		let yhat = predict(&imodel, &predictstate, errarray, &predictvars);
		
		printf("%f %f %f %f %f\n", ptr->ID, ptr->TIME, ptr->DV, yhat, predictvars.IPRED);
	}
	advanfuncs->destruct(advan);
	free(advan);
	advanfuncs_free(advanfuncs);

	return EXIT_SUCCESS;
}

