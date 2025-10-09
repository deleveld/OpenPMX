/*
 * Theophylline oral from Bae and Yim from TCP 2016
 * compile and run with:
 * gcc -W -Wall -Wextra -O2 theo_posthoc.c -I../../include -I../../src -lgsl -lgslcblas -lm; ./a.out
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

extern RECORD data[11]; /* forward declaration */

/* include sources directly instead of linking to library */
#include "dataconfig/dataconfig.c"
#include "advan/advan.c"
#include "advan/pred.c"
#include "bobyqa/bobyqa.c"
#include "idata.c"
#include "print.c"
#include "options.c"
#include "popmodel.c"
#include "pmxstate.c"
#include "omegainfo.c"
#include "linalg.c"
#include "ievaluate.c"
#include "predict.c"
#include "stage1.c"
#include "scatter.c"
#include "utils/vector.c"
#include "utils/various.c"

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

	let advanfuncs = advanfuncs_alloc(&openpmx.data, &openpmx.advan);
	let popmodel = popmodel_init(&openpmx);
	var idata = idata_construct(&advanfuncs->recordinfo,
								popmodel.ntheta,
								popmodel.nomega,
								popmodel.nsigma,
								advanfuncs->nstate,
								advanfuncs->advanconfig->imodelfields.size,
								advanfuncs->advanconfig->predictfields.size);
	let omegainfo = omegainfo_init(popmodel.nomega,
								   popmodel.omega,
								   popmodel.omegafixed);
	let o = (OPTIONS){ };
	let options = options_default(&o);

	stage1_thread(idata.individ,
				  advanfuncs,
				  &popmodel,
				  &omegainfo.nonzero,
				  &options,
				  0 /*scatteroptions*/);

	idata_predict_pred_thread(idata.individ,
				  advanfuncs,
				  &popmodel,
				  0,
				  &options,
				  0 /*scatteroptions*/);

	forcount(i, idata.nindivid) {
		let individ = &idata.individ[i];
		let imodel = &individ->imodel[i];

		printf("ID %f\n", individ->ID);
		printf("V %f\n", imodel->V);
		printf("K %f\n", imodel->K);
		printf("KA %f\n", imodel->KA);
	
		forcount(n, individ->nrecord) {
			let record = &individ->record[n];
			let ipred = individ->yhat[n];
			let pred = individ->pred[n];
			printf("%f %f %f %f %f\n", record->ID, record->TIME, record->DV, ipred, pred);
		}
	}

	advanfuncs_free(advanfuncs);
	idata_destruct(&idata);
	
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

