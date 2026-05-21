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

/// This file implements the Eleveld remifentanil model.
///
/// Eleveld DJ, Proost JH, Vereecke H, Absalom AR, Olofsen E, Vuyk J,
/// Struys MM. An allometric model of remifentanil pharmacokinetics and
/// pharmacodynamics. Anesthesiology. 2017 Jun 1;126(6):1005-18.
///

#include <math.h>

#include "openpmx_model.h"

static double remifentanil_eleveld_theta[] = {
	 1.759110e+00, //  v1=5.81
	 2.177170e+00, //  v2=8.82
	 1.614450e+00, //  v3=5.03
	 9.462570e-01, //  cl=2.58
	 5.403150e-01, //  q2=1.72
	-2.083880e+00, //  q3=0.12
	 2.878980e+00, //  e50 for CL maturation
	-5.544810e-03, //  aging v1/q2/q3
	-3.269850e-03, //  aging v2/cl
	-3.151350e-02, //  aging v3
	 4.704050e-01, //  increase cl/q2/v2 in females 12-45 years
	-2.604960e-02, //  weight correction v3
	0, // (0.01, 1.113080e-01, 1) ; residual error Minto data
	0, // (0.01, 2.714780e-01, 1) ; residual error Ross data
	0, // (0.01, 2.402460e-01, 1) ; residual error Mertens data
};

static double THETA(const int i)
{
	return remifentanil_eleveld_theta[i - 1];
}

static double ETA(const int i)
{
	(void) i;
	return 0.;
}

REMIFENTANIL_ELEVELD pmx_model_remifentanil_eleveld(REMIFENTANIL_ELEVELD_COVARIATES* config)
{
	const double AGE = config->age;
	const double WGT = config->weight;
	const double HGT = config->height;
	const double M1F2 = config->female ? 2 : 1;

	// maturation
	const double SE50=THETA(7);
	const double ADLT=pow(WGT,2)/(pow(WGT,2)+pow(SE50,2));
	const double AREF=pow(70.,2)/(pow(70.,2)+pow(SE50,2));
	const double KMAT=ADLT/AREF;
	// scaling using Al-sallami FFM
	const double HT2=(HGT/100.)*(HGT/100.);
	const double MATM=0.88+((1-0.88)/(1+pow(AGE/13.4,-12.7)));
	const double MATF=1.11+((1-1.11)/(1+pow(AGE/7.1,-1.1)));
	const double FFMM=MATM*42.92*(HT2)*WGT/(30.93*(HT2)+WGT);
	const double FFMF=MATF*37.99*(HT2)*WGT/(35.98*(HT2)+WGT);
	const double FFMR=42.92*(1.7*1.7)*70./(30.93*(1.7*1.7)+70.);
	const double MAL=2-M1F2;
	const double FEM=M1F2-1;
	const double BSIZ=(MAL*FFMM + FEM*FFMF)/FFMR;
	// aging for v1/q2/q2, v3 and v2/cl
	const double KV1=exp(THETA(8)*(AGE-35.));
	const double KV2=exp(THETA(9)*(AGE-35.));
	const double KV3=exp(THETA(10)*(AGE-35.));
	const double KCL=KV2;
	const double KQ2=KV1;
	const double KQ3=KV1;
	// sex correction for cl, v2 and q2
	const double PPUB=(pow(AGE,6))/(pow(AGE,6) + pow(12,6));
	const double ELDY=(pow(AGE,6))/(pow(AGE,6) + pow(45,6));
	const double KSEX=1+(M1F2-1)*PPUB*(1-ELDY)*THETA(11);
	// weight correction for V3;
	const double WV3=exp(THETA(12)*(WGT-70.));
	// seperate residual error for each study
//	const double RESV=0;
	// compartmental allometric scaling
	const double M1 =BSIZ * KV1;
	const double M2 =BSIZ * KV2 * KSEX;
	const double M3 =BSIZ * KV3 * WV3;
	const double V1 =exp(THETA(1)+ETA(1)) * M1;
	const double V2 =exp(THETA(2)+ETA(2)) * M2;
	const double V3 =exp(THETA(3)+ETA(3)) * M3;
	const double RV2=exp(THETA(2));
	const double RV3=exp(THETA(3));
	const double M4 =pow(BSIZ,0.75) * KCL * KSEX * KMAT;
	const double M5 =pow(V2/RV2,0.75) * KQ2 * KSEX;
	const double M6 =pow(V3/RV3,0.75) * KQ3;
	const double CL =exp(THETA(4)+ETA(4)) * M4;
	const double Q2 =exp(THETA(5)+ETA(5)) * M5;
	const double Q3 =exp(THETA(6)+ETA(6)) * M6;
	const double KE0=1.09*exp(-0.0289*(AGE-35));

	return (REMIFENTANIL_ELEVELD) {
		.FFM = MAL*FFMM + FEM*FFMF,
		.V1 = V1,
		.V2 = V2,
		.V3 = V3,
		.CL = CL,
		.Q2 = Q2,
		.Q3 = Q3,
		.KE0 = KE0,
		.k10 = CL / V1,
		.k12 = Q2 / V1,
		.k21 = Q2 / V2,
		.k13 = Q3 / V1,
		.k31 = Q3 / V3,
	};
}
