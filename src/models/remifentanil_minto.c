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

/// This file implements the Minto remifentanil model.
///
/// Minto CF, Schnider TW, Egan TD, Youngs E, Lemmens HJ, Gambus PL,
/// Billard V, Hoke JF, Moore KH, Hermann DJ, Muir KT. Influence of age
/// and gender on the pharmacokinetics and pharmacodynamics of
/// remifentanil: I. Model development. Anesthesiology. 1997 Jan 1;86(1):10-23.
///
/// And:
/// Minto CF, Schnider TW, Shafer SL. Pharmacokinetics and pharmacodynamics
/// of remifentanil: II. Model application. Anesthesiology. 1997 Jan
/// 1;86(1):24-33.
///

#include "openpmx_model.h"

REMIFENTANIL_MINTO pmx_model_remifentanil_minto(REMIFENTANIL_MINTO_COVARIATES* config)
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

	const double V1 = 5.1 - 0.0201 * (age - 40.) + 0.072 * (LBM - 55.);
	const double V2 = 9.82 - 0.0811 * (age - 40.) + 0.108 * (LBM - 55.);
	const double V3 = 5.42;
	const double CL = 2.6 - 0.0162 * (age - 40.) + 0.0191 * (LBM - 55.);
	const double Q2 = 2.05 - 0.0301 * (age - 40.);
	const double Q3 = 0.076 - 0.00113 * (age - 40.);
	const double KE0 = 0.595 - 0.007 * (age - 40.);

	return (REMIFENTANIL_MINTO) {
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
