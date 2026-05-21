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

#ifndef OPENPMX_MODELH
#define OPENPMX_MODELH

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------------------------------------------------------*/
typedef struct {
	const double age;		/* years */
	const double weight;	/* kg */
	const double height;	/* cm */
	const bool female;
} PROPOFOL_SCHNIDER_COVARIATES;

typedef struct {
	const double LBM;
	const double V1, V2, V3, CL, Q2, Q3, KE0;
	const double k10, k12, k21, k13, k31;
} PROPOFOL_SCHNIDER;

PROPOFOL_SCHNIDER pmx_model_propofol_schnider(PROPOFOL_SCHNIDER_COVARIATES* config);

/*---------------------------------------------------------------------*/
typedef struct {
	const double age;		/* years */
	const double weight;	/* kg */
	const double height;	/* cm */
	const bool female;
	const bool opiates;
} PROPOFOL_ELEVELD_COVARIATES;

typedef struct {
	const double FFM;
	const double E50;
	const double EMAX;
	const double GAM;
	const double GAM1;
	const double V1, V2, V3, CL, Q2, Q3, KE0;
	const double k10, k12, k21, k13, k31;
} PROPOFOL_ELEVELD;

PROPOFOL_ELEVELD pmx_model_propofol_eleveld(PROPOFOL_ELEVELD_COVARIATES* config);

double pmx_model_propofol_eleveld_bis(PROPOFOL_ELEVELD* model, const double ceff);

/*---------------------------------------------------------------------*/
typedef struct {
	const double age;		/* years */
	const double weight;	/* kg */
	const double height;	/* cm */
	const bool female;
} REMIFENTANIL_ELEVELD_COVARIATES;

typedef struct {
	const double FFM;
	const double V1, V2, V3, CL, Q2, Q3, KE0;
	const double k10, k12, k21, k13, k31;
} REMIFENTANIL_ELEVELD;

REMIFENTANIL_ELEVELD pmx_model_remifentanil_eleveld(REMIFENTANIL_ELEVELD_COVARIATES* config);

/*---------------------------------------------------------------------*/
typedef struct {
	const double age;		/* years */
	const double weight;	/* kg */
	const double height;	/* cm */
	const bool female;
} REMIFENTANIL_MINTO_COVARIATES;

typedef struct {
	const double LBM;
	const double V1, V2, V3, CL, Q2, Q3, KE0;
	const double k10, k12, k21, k13, k31;
} REMIFENTANIL_MINTO;

REMIFENTANIL_MINTO pmx_model_remifentanil_minto(REMIFENTANIL_MINTO_COVARIATES* config);

/*---------------------------------------------------------------------*/
typedef struct {
	const double V1, V2, V3, CL, Q2, Q3;
	const double k10, k12, k21, k13, k31;
} SUFENTANIL_GEPTS;

SUFENTANIL_GEPTS pmx_model_sufentanil_gepts(void);

/*---------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif
