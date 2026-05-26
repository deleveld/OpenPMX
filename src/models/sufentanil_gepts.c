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

/// This file implements the Gepts sufentanil model
///
/// Gepts E, Shafer SL, Camu F, Stanski DR, Woestenborghs R, Van Peer A,
/// Heykants JJ. Linearity of pharmacokinetics and model estimation of
/// sufentanil. Anesthesiology. 1995 Dec 1;83(6):1194-204. 
///

#include "openpmx_model.h"

SUFENTANIL_GEPTS pmx_model_sufentanil_gepts(void)
{
	const double V1 = 14.3;
	const double V2 = 63.1;
	const double V3 = 261.6;
	const double CL = 0.92;
	const double Q2 = 1.55;
	const double Q3 = 0.33;
	
	return (SUFENTANIL_GEPTS) {
		.V1 = V1,
		.V2 = V2,
		.V3 = V3,
		.CL = CL,
		.Q2 = Q2,
		.Q3 = Q3,
		.k10 = CL / V1,
		.k12 = Q2 / V1,
		.k21 = Q2 / V2,
		.k13 = Q3 / V1,
		.k31 = Q3 / V3,
	};
}
