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

#include <math.h>

#include "openpmx.h"

static double eleveld_theta[] = {
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
};

#define THETA(i)	eleveld_theta[(i)-1]
#define ETA(i) 		0.

TCIMODEL pmx_advan_tci_eleveld_propofol(const double AGE,
										const double WGT,
										const double HGT,
										const bool female,
										const bool opiates)
{
	const double A1V2 = 1.; /* assume arterial */
	const double M1F2 = (female) ? 2. : 1.;
	const double PMA = AGE + 40./52;
	const double TECH = (opiates) ? 2. : 1.;

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

	const double TKE0 = -1.921620e+00;
	const double ke0 = exp(TKE0+ETA(12))*pow(WGT/70.,-0.25);

	const double k10 = CL / V1;
	const double k12 = Q2 / V1;
	const double k21 = Q2 / V2;
	const double k13 = Q3 / V1;
	const double k31 = Q3 / V3;

	return (TCIMODEL) {
		.k10 = k10,
		.k12 = k12,
		.k13 = k13,
		.k21 = k21,
		.k31 = k31,
		.ke0 = ke0,
		.vc = V1,
	};
}
