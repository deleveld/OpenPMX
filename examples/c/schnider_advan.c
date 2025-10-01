/*
 * Prediction from single individual with a three compartment model
 * 
 * PK model and data from first individual:
 * Schnider T, Minto C, Gambus P, Andresen C, Goodale D, Shafer S, Youngs E:
 * The influence of method of administration and covariates on the pharmacokinetics
 * of propofol in adult volunteers. Anesthesiology 1998; 88:1170–82 PMID: 9605675
 * 
 * Compile and run with:
 * gcc -W -Wall -Wextra -O2 schnider_advan.c -I../../include -I../../src -lgsl -lgslcblas -lm; ./a.out
 */

#include <stdlib.h>
#include <math.h>
#include <stddef.h>

#include "openpmx.h"

typedef struct IMODEL {
	double V1;
	double V2;
	double V3;
	double CL;
	double Q2;
	double Q3;
} IMODEL;

typedef struct RECORD {
	double ID, TIME, DV, AMT, RATE, EVID, AGE, WT, HT, M1F2, A1V2;
} RECORD;

static void imodel_init(IMODEL* const _imodel,
						ADVANSTATE* const _advanstate,
						const POPPARAM* const _popparam)
{
	const RECORD* const _record = _advanstate->record;
	const double WGT = _record->WT;
	const double AGE = _record->AGE;
	const double HGT = _record->HT;
	const double M1F2 = _record->M1F2;
	(void)_popparam;

/* james equation */
    const double LBMM=1.1 *WGT-128*(WGT/HGT)*(WGT/HGT);
    const double LBMF=1.07*WGT-148*(WGT/HGT)*(WGT/HGT);
    const double FEM =M1F2-1;
    const double LBM =(1-FEM)*LBMM+(FEM)*LBMF;
/* schnider model */
    _imodel->V1 = 4.27;
    _imodel->V2 = 18.9-0.391*(AGE-53);
    _imodel->V3 = 238;
    _imodel->CL = 1.89+0.0456*(WGT-77)-0.0681*(LBM-59)+0.0264*(HGT-177);
    if (_imodel->CL <= 0.) 
		_imodel->CL = 0.01;
    _imodel->Q2 = 1.29-0.024*(AGE-53);
    _imodel->Q3 = 0.836;
}

typedef struct PREDICTVARS {
	double IPRED;
} PREDICTVARS;

static void imodel_diffeqn(double _dadt[],
						   const IMODEL* const _imodel,
						   const RECORD* const _record,
						   const double* const _state,
						   const POPPARAM* const _popparam,
						   const double T)
{
	(void) _record;
	(void) _popparam;
	(void) T;
	
	const double V1 = _imodel->V1;
	const double V2 = _imodel->V2;
	const double V3 = _imodel->V3;
	const double CL = _imodel->CL;
	const double Q2 = _imodel->Q2;
	const double Q3 = _imodel->Q3;

#define THETA(i) 	((const double)_theta[i-1])
#define ETA(i) 		((const double)_eta[i-1])
#define A(i) 		((const double)_state[i-1])
#define DADT(i) 	(_dadt[i-1])
	const double K10 = CL/V1;
	const double K12 = Q2/V1;
	const double K21 = Q2/V2;
	const double K13 = Q3/V1;
	const double K31 = Q3/V3;
	const double K1023 = K10 + K12 + K13;
	DADT(1) = A(2)*K21 + A(3)*K31 - A(1)*K1023;
	DADT(2) = A(1)*K12 - A(2)*K21;
	DADT(3) = A(1)*K13 - A(3)*K31;
#undef THETA
#undef ETA
#undef A
#undef DADT
}

static double imodel_predict(const IMODEL* const _imodel,
							 const PREDICTSTATE* const _predictstate,
							 const POPPARAM* const _popparam,
							 const double* const _err,
							 PREDICTVARS* _predparams)
{
	(void)_popparam;
	
	const double* const _state = _predictstate->state;
	double Y = NAN;
	const double V1 = _imodel->V1;
	double IPRED;

#define A(i) 	((const double)_state[i-1])
#define ERR(i) 	((const double)_err[i-1])
	IPRED = A(1)/V1;
	Y = IPRED*(1 + ERR(1));
#undef A
#undef ERR

	_predparams->IPRED = IPRED;

	return Y;
}

#define ARRAYSIZE(a) (sizeof(a)/sizeof(a[0]))

/* ID, TIME, DV, AMT, RATE, EVID, AGE, WT, HT, M1F2, A1V2 */
const RECORD data[] = {
	{	1,0,0,92.60001,342.963,4,34,46.3,157.5,2,1	},
	{	1,2.11,3.62,0,0,0,34,46.3,157.5,2,1	},
	{	1,4.01,1.33,0,0,0,34,46.3,157.5,2,1	},
	{	1,8.06,0.676,0,0,0,34,46.3,157.5,2,1	},
	{	1,16.03,0.338,0,0,0,34,46.3,157.5,2,1	},
	{	1,30.05,0.204,0,0,0,34,46.3,157.5,2,1	},
	{	1,59.73,0.11,0,0,0,34,46.3,157.5,2,1	},
	{	1,60,0,89.4,1.49,1,34,46.3,157.5,2,1	},
	{	1,62,0.322,0,0,0,34,46.3,157.5,2,1	},
	{	1,63.96,0.632,0,0,0,34,46.3,157.5,2,1	},
	{	1,67.99,0.561,0,0,0,34,46.3,157.5,2,1	},
	{	1,76.18,0.606,0,0,0,34,46.3,157.5,2,1	},
	{	1,90.02,0.748,0,0,0,34,46.3,157.5,2,1	},
	{	1,119.68,0.812,0,0,0,34,46.3,157.5,2,1	},
	{	1,122.03,0.373,0,0,0,34,46.3,157.5,2,1	},
	{	1,124.11,0.291,0,0,0,34,46.3,157.5,2,1	},
	{	1,128.06,0.23,0,0,0,34,46.3,157.5,2,1	},
	{	1,136,0.184,0,0,0,34,46.3,157.5,2,1	},
	{	1,150.03,0.16,0,0,0,34,46.3,157.5,2,1	},
	{	1,180,0.111,0,0,0,34,46.3,157.5,2,1	},
	{	1,240,0.0864,0,0,0,34,46.3,157.5,2,1	},
	{	1,300,0.0509,0,0,0,34,46.3,157.5,2,1	},
	{	1,598,0.027,0,0,0,34,46.3,157.5,2,1	},
	{	1,0,0,92.60001,342.963,4,34,46.3,157.5,2,1	},
	{	1,2.01,2.63,0,0,0,34,46.3,157.5,2,1	},
	{	1,4,1.17,0,0,0,34,46.3,157.5,2,1	},
	{	1,8.01,0.567,0,0,0,34,46.3,157.5,2,1	},
	{	1,16.07,0.246,0,0,0,34,46.3,157.5,2,1	},
	{	1,30,0.151,0,0,0,34,46.3,157.5,2,1	},
	{	1,59.58,0.127,0,0,0,34,46.3,157.5,2,1	},
	{	1,60,0,89.4,1.49,1,34,46.3,157.5,2,1	},
	{	1,61.98,0.379,0,0,0,34,46.3,157.5,2,1	},
	{	1,64.01,0.433,0,0,0,34,46.3,157.5,2,1	},
	{	1,68.01,0.506,0,0,0,34,46.3,157.5,2,1	},
	{	1,76.02,0.62,0,0,0,34,46.3,157.5,2,1	},
	{	1,90,0.501,0,0,0,34,46.3,157.5,2,1	},
	{	1,119.55,0.609,0,0,0,34,46.3,157.5,2,1	},
	{	1,121.96,0.294,0,0,0,34,46.3,157.5,2,1	},
	{	1,124.01,0.233,0,0,0,34,46.3,157.5,2,1	},
	{	1,128.02,0.191,0,0,0,34,46.3,157.5,2,1	},
	{	1,136,0.152,0,0,0,34,46.3,157.5,2,1	},
	{	1,149.98,0.121,0,0,0,34,46.3,157.5,2,1	},
	{	1,180.2,0.1,0,0,0,34,46.3,157.5,2,1	},
	{	1,240,0.051,0,0,0,34,46.3,157.5,2,1	},
	{	1,300,0.044,0,0,0,34,46.3,157.5,2,1	},
	{	1,600.47,0.026,0,0,0,34,46.3,157.5,2,1	},
};

/* include sources directly instead of linking to library */
#include "dataconfig/dataconfig.c"
#include "advan/advan.c"
#include "advan/diffeqn_libgsl.c"
#include "print.c"
#include "popmodel.c"
#include "utils/vector.c"

int main(void)
{
	OPENPMX openpmx = (OPENPMX) {
		.data = {
			.records = data,
			.nrecords = ARRAYSIZE(data),
			.recordfields = {
				.size = sizeof(RECORD),
				.field = {
					{ .name="ID", .offset=offsetof(RECORD, ID) },
					{ .name="TIME", .offset=offsetof(RECORD, TIME) },
					{ .name="AMT", .offset=offsetof(RECORD, AMT) },
					{ .name="RATE", .offset=offsetof(RECORD, RATE) },
					{ .name="EVID", .offset=offsetof(RECORD, EVID) },
					{ .name="WT", .offset=offsetof(RECORD, WT) },
				},
			},
		},
		.advan = {
			.init = imodel_init,
			.predict = imodel_predict,
			.diffeqn = imodel_diffeqn,
			.method = pmx_advan_diffeqn_libgsl,
			.firstonly = 1,
			.nstate = 3,
		},
		.theta = { },
		.omega = { },
		.sigma = { },
	};

	/* get advancer function table, memory needed to iterate, and construct advancer */
	let advanfuncs = advanfuncs_alloc(&openpmx.data, &openpmx.advan);
	ADVAN* advan = calloc(advanfuncs->advan_size, 1);	
	advanfuncs->construct(advan, advanfuncs);
	
	/* arguments needed to iterate */
	double eta[OPENPMX_OMEGA_MAX] = { };
	let popmodel = popmodel_init(&openpmx);
	let popparam = (POPPARAM) { 
		.theta = popmodel.theta,
		.ntheta = popmodel.ntheta,
		.eta = eta,
		.nomega = popmodel.nomega,
		.sigma = popmodel.sigma,
		.nsigma = popmodel.nsigma,
		.nstate = advanfuncs->nstate,
	};
	
	/* scratch memory needed to iterate */
	double errarray[OPENPMX_SIGMA_MAX] = { };
	IMODEL imodel;
	PREDICTVARS predictvars;

	/* advance over individuals data */
	let predict = openpmx.advan.predict;
	forcount(i, ARRAYSIZE(data)) {
		let ptr = &data[i];
		let predictstate = advan_advance(advan, &imodel, ptr, &popparam);
		let yhat = predict(&imodel, &predictstate, &popparam, errarray, &predictvars);
		
		printf("%f %f %f %f %f A1=%f A2=%f A3=%f\n", ptr->ID, ptr->TIME, ptr->DV, yhat, predictvars.IPRED,
													predictstate.state[0],
													predictstate.state[1],
													predictstate.state[2]);
	}
	advanfuncs->destruct(advan);
	free(advan);

	return EXIT_SUCCESS;
}

