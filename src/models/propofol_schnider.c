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

/// This file implements the Schnider propofol model.
///
///	Schnider T, Minto C, Gambus P, Andresen C, Goodale D, Shafer S, Youngs E:
///	The influence of method of administration and covariates on the pharmacokinetics
///	of propofol in adult volunteers. Anesthesiology 1998; 88:1170–82 PMID: 9605675
///

#include "openpmx_model.h"

PROPOFOL_SCHNIDER pmx_model_propofol_schnider(PROPOFOL_SCHNIDER_COVARIATES* config)
{
	const double age = config->age;
	const double weight = config->weight;
	const double height = config->height;
	const double m1f2 = config->female ? 2 : 1;

	double LBM = 0.;
	if (m1f2 == 1)
		LBM = 1.1 * weight - 128. * (weight/height) * (weight/height);
	else
		LBM = 1.07 * weight - 148. * (weight/height) * (weight/height);

	const double V1 = 4.27;
	const double V2 = 18.9 - 0.391 * (age - 53.);
	const double V3 = 238.;
	const double CL = 1.89 + (weight - 77) * 0.0456 - (LBM - 59) * 0.0681 + (height - 177.) * 0.0264;
	const double Q2 = 1.29 - 0.024 * (age - 53.);
	const double Q3 = 0.836;
	const double KE0 = 0.456;

	return (PROPOFOL_SCHNIDER) {
		.LBM = LBM,
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
