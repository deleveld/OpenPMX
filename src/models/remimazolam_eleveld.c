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

/// This file implements the Eleveld remimazolam model.
///
/// Eleveld DJ, Colin PJ, Van den Berg JP, Koomen JV, Stoehr T, Struys MM. 
/// Development and analysis of a remimazolam pharmacokinetics and 
/// pharmacodynamics model with proposed dosing and concentrations for 
/// anaesthesia and sedation. British Journal of Anaesthesia. 
/// 2025 Jul 1;135(1):206-17.
///

#include "openpmx_model.h"

#include <math.h>

REMIMAZOLAM_ELEVELD pmx_model_remimazolam_eleveld(REMIMAZOLAM_ELEVELD_COVARIATES* config)
{
	const double AGE = config->age;
	const double WGT = config->weight;
	const double M1F2 = config->female ? 2 : 1;
	const double OA1P2 = config->opiates ? 2. : 1.;

	const double KOCL = -0.139340*(OA1P2-1.);
	const double KSV3 = 0.287037*(M1F2-1.);
	const double KSCL = 0.162802*(M1F2-1.);
	const double KAV3 = 0.00730817*(AGE-35.);

	const double VSIZ = WGT / 70.;
	const double CSIZ = pow(WGT / 70., 0.75);
	const double V1 = VSIZ * 4.30730;
	const double V2 = VSIZ * 12.2994;
	const double V3 = VSIZ * 18.6411 * exp(KAV3 + KSV3);
	const double CL = CSIZ * 1.11977 * exp(KSCL + KOCL);
	const double Q2 = pow(VSIZ, 0.75) * 1.45260;
	const double Q3 = pow(V3 / 18.6411, 0.75) * 0.297838;
	const double KE0 = 0.298269;
	const double E50 = 0.182 * exp(-7.63/1000. * (AGE - 35.));

	return (REMIMAZOLAM_ELEVELD) {
		.V1 = V1,
		.V2 = V2,
		.V3 = V3,
		.CL = CL,
		.Q2 = Q2,
		.Q3 = Q3,
		.KE0 = KE0,
		.E50 = E50,
		.k10 = CL / V1,
		.k12 = Q2 / V1,
		.k21 = Q2 / V2,
		.k13 = Q3 / V1,
		.k31 = Q3 / V3,
	};
}
