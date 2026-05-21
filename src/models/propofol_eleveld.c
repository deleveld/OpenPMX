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

/// This file implements the Eleveld propofol model.
///
/// Eleveld DJ, Colin P, Absalom AR, Struys MM. Pharmacokinetic–pharmacodynamic
/// model for propofol for broad application in anaesthesia and sedation.
/// British journal of anaesthesia. 2018 May 1;120(5):942-59.
///

#include <math.h>

#include "openpmx_model.h"

static double propofol_eleveld_theta[] = {
	1.837860e+00,
	3.238730e+00,
	5.608800e+00,
	5.819830e-01,
	5.596720e-01,
	1.030460e-01,
	1.913070e-01,
	3.744220e+00,
	2.203300e+00,
	-1.563300e-02,
	-2.857090e-03,
	3.513130e+00,
	-1.381660e-02,
	4.223570e+00,
	7.420430e-01,
	2.656420e-01,
	3.498850e-01,
	-3.849270e-01,
	0,
	0,
	1.124750e+00,	// e50=3.08 mg/l
	-1.921620e+00,	// ke0=0.146 1/min
	9.298240e+01,	// emax
	3.877210e-01,	// gamma=1.47
	2.082830e+00,	// residual error BIS=8.03
	5.173990e-02,	// age delay
	-6.348370e-03,	// age e50
	2.156050e-01,	// venous ke0=1.24 1/min
	6.386410e-01,	// gamma=1.89
};

static double THETA(const int i)
{
	return propofol_eleveld_theta[i - 1];
}

static double ETA(const int i)
{
	(void) i;
	return 0.;
}

PROPOFOL_ELEVELD pmx_model_propofol_eleveld(PROPOFOL_ELEVELD_COVARIATES* config)
{
	const double AGE = config->age;
	const double WGT = config->weight;
	const double HGT = config->height;
	const double M1F2 = config->female ? 2 : 1;
	
	const double A1V2 = 1.; 			/* assume arterial */
	const double PMA = AGE + 40./52;	/* assume full-term */
	const double TECH = config->opiates ? 2. : 1.;

	// Al-sallami FFM
	const double HT2=(HGT/100.)*(HGT/100.);
	const double MATM=0.88+((1-0.88)/(1+pow(AGE/13.4,-12.7)));
	const double MATF=1.11+((1-1.11)/(1+pow(AGE/7.1,-1.1)));
	const double MATR=0.88+((1-0.88)/(1+pow(35./13.4,-12.7)));
	const double FFMM=MATM*42.92*(HT2)*WGT/(30.93*(HT2)+WGT);
	const double FFMF=MATF*37.99*(HT2)*WGT/(35.98*(HT2)+WGT);
	const double FFMR=MATR*42.92*(1.7*1.7)*70./(30.93*(1.7*1.7)+70.);
	const double FFM=FFMM*(2-M1F2) + FFMF*(M1F2-1);
	const double NFFM=FFM/FFMR;
	// maturation
	const double DV1=1;
	const double DV2=1;
	const double DV3=1;
	// sigmoidal maturation of CL
	const double PMW=PMA*52.;
	const double PMR=(35.+40./52.)*52.;
	const double ME50=exp(THETA(8));
	const double MGAM=exp(THETA(9));
	const double MCL=pow(PMW,MGAM)/(pow(PMW,MGAM)+pow(ME50,MGAM));
	const double RCL=pow(PMR,MGAM)/(pow(PMR,MGAM)+pow(ME50,MGAM));
	const double DCL=MCL/RCL;
	const double DQ2=1;
	// sigmoidal maturation of Q3 based on 40 weeks gestation
	const double PMEW=AGE*52.+40.;
	const double PMER=35.*52.+40.;
	const double QE50=exp(THETA(14));
	const double MQ3=PMEW/(PMEW+QE50);
	const double RQ3=PMER/(PMER+QE50);
	const double DQ3=MQ3/RQ3;
	// aging
	const double KV1=1;
	const double KV2=exp(THETA(10)*(AGE-35.));
	const double KV3=exp(THETA(13)*(AGE)*(TECH-1));
	const double KCL=exp(THETA(11)*(AGE)*(TECH-1));
	const double KQ2=1;
	const double KQ3=1;
	// covariate structure
	// V1 scales sigmoid with weight
	const double VV50=exp(THETA(12));
	const double CV1=WGT/(WGT+VV50);
	const double RV1=70./(70.+VV50);
	const double M1 =(CV1/RV1) * KV1 * DV1;
	const double VCV1=(A1V2-1)*(1-CV1);
	const double V1 =exp(THETA(1)+ETA(1)) * M1 * (1+VCV1*exp(THETA(17)));
	const double M2 =(WGT/70.) * KV2 * DV2;
	const double V2 =exp(THETA(2)+ETA(2)) * M2;
	const double M3 =(NFFM) * KV3 * DV3;
	const double V3 =exp(THETA(3)+ETA(3)) * M3;
	const double M4 =pow(WGT/70.,0.75) * KCL * DCL;
	const double CL =exp((2-M1F2)*THETA(4)+(M1F2-1)*THETA(15)+ETA(4)) * M4;
	const double RV2=exp(THETA(2));
	const double M5 =pow(V2/RV2,0.75) * KQ2 * DQ2;
	const double KM5=1+exp(THETA(16))*(1-MQ3);
	const double Q2 =exp(THETA(5)+ETA(5)+(A1V2-1)*THETA(18)) * M5 * KM5;
	const double RV3=exp(THETA(3));
	const double M6 =pow(V3/RV3,0.75) * KQ3 * DQ3;
	const double Q3 =exp(THETA(6)+ETA(6)) * M6;

	// effect compartment
	const double E50 =exp(THETA(21)+ETA(11)+THETA(27)*(AGE-35.));
	const double TKE0=(2-A1V2)*THETA(22) + (A1V2-1)*THETA(28);
	const double KE0 =exp(TKE0+ETA(12))*pow(WGT/70.,-0.25);
	const double EMAX=THETA(23);
	const double GAM =exp(THETA(24));
	const double GAM1=exp(THETA(29));

	return (PROPOFOL_ELEVELD) {
		.FFM = FFM,
		.V1 = V1,
		.V2 = V2,
		.V3 = V3,
		.CL = CL,
		.Q2 = Q2,
		.Q3 = Q3,
		.KE0 = KE0,
		.E50 = E50,
		.EMAX = EMAX,
		.GAM = GAM,
		.GAM1 = GAM1,
		.k10 = CL / V1,
		.k12 = Q2 / V1,
		.k21 = Q2 / V2,
		.k13 = Q3 / V1,
		.k31 = Q3 / V3,
	};
}

double pmx_model_propofol_eleveld_bis(PROPOFOL_ELEVELD* model, const double ceff)
{
	const double e50 = model->E50;
	const double emax = model->EMAX;
	const double gam = model->GAM;
	const double gam1 = model->GAM1;

	/* low concentrations (< E50) get steeper gamma */
	const double sig_low  = pow(ceff, gam1) / (pow(ceff, gam1) + pow(e50, gam1));
	const double sig_high = pow(ceff, gam) / (pow(ceff, gam) + pow(e50, gam));
	if (ceff < e50)
		return emax * (1. - sig_low);
	return emax * (1. - sig_high);
}

