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

#ifndef OPENPMX_H
#define OPENPMX_H

#include <stdbool.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------------------------------------------------------*/
/// This defines the version number and the various maximum size of
/// theta, omega, and sigma. Also the maximum number of fields in a record,
/// the number of state variables etc. If these are too low they can
/// be adjusted here and the library must be recompiled.
#define OPENPMX_VERSION_MAJOR			0
#define OPENPMX_VERSION_MINOR			1
#define OPENPMX_VERSION_RELEASE			2

#define OPENPMX_THETA_MAX				64
#define OPENPMX_OMEGABLOCKSIZE_MAX		64
#define OPENPMX_OMEGABLOCK_MAX			10
#define OPENPMX_OMEGA_MAX				16
#define OPENPMX_SIGMA_MAX				8

#define OPENPMX_FIELDS_MAX				64
#define OPENPMX_FIELDNAME_MAX			64

/*---------------------------------------------------------------------*/
/* dataconfig */
/*---------------------------------------------------------------------*/
typedef	struct {
	int size;
	struct {
		char name[OPENPMX_FIELDNAME_MAX];
		int offset;
	} field[OPENPMX_FIELDS_MAX];
} STRUCTINFO;

typedef struct RECORD RECORD;
typedef struct {
	RECORD* writeable;
	const RECORD* records;
	const int nrecords;
	const bool _offset1;
	const STRUCTINFO recordfields;
} DATACONFIG;

/*---------------------------------------------------------------------*/
/* advan */
/*---------------------------------------------------------------------*/
typedef	struct {
	const double* const theta;
	const int ntheta;
	const double* const eta;
	const int nomega;
	const double* const sigma;
	const int nsigma;
	const int nstate;
} POPPARAM;

/* info about the record being processed */
typedef struct {
	const double statetime;
	const double* const state;
	const RECORD* const record;
	const POPPARAM* const popparam;
} PREDICTSTATE;

/* info about the state of the advancer, only available in the init */
typedef struct ADVAN ADVAN;
typedef struct {
	ADVAN* const advan;
	const PREDICTSTATE current;
} ADVANSTATE;

/* changes to advancer within the IMODEL_INIT function */
void pmx_advan_amtlag(const ADVANSTATE* advanstate, const int cmt, const double t);
void pmx_advan_bioaval(const ADVANSTATE* advanstate, const int cmt, const double f);
bool pmx_advan_inittime(const ADVANSTATE* advanstate, const double t);
void pmx_advan_state_init(const ADVANSTATE* advanstate, const int cmt, const double v);
void pmx_advan_eigen_sysmat(const ADVANSTATE* advanstate, const double* sysmat);

/* callback for differential equation solver */
typedef struct IMODEL IMODEL;
typedef struct PREDICTVARS PREDICTVARS;
typedef struct ADVANCONFIG ADVANCONFIG;
typedef struct ADVANFUNCS ADVANFUNCS;
typedef struct ADVANCONFIG {
	
	/* init function */
	void (*init)(IMODEL* const imodel, 
				 ADVANSTATE* advanstate);
	/* predict function */
	double (*predict)(const IMODEL* const imodel,
					  const PREDICTSTATE* const predictstate,
					  const double* const err,
					  PREDICTVARS* predparams);
	const int firstonly;
	const int predictall;

	const STRUCTINFO imodelfields;
	const STRUCTINFO predictfields;
	ADVANFUNCS* (*method)(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);
	const int nstate;
	
	/* differential equation callback */
	void (*diffeqn)(double DADT[],
					const IMODEL* const imodel,
					const RECORD* const record,
					const double* state,
					const POPPARAM* const popparam,
					const double T);

	/* extra arguments */
	union {
		/* for diffeqn solvers */
		struct {
			const double hstart;
			const double abstol;
			const double reltol;
			const char* steptype;
		} diffeqn;
	} args;
} ADVANCONFIG;

/* basic analytic advan methods */
ADVANFUNCS* pmx_advan_pred(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);
ADVANFUNCS* pmx_advan_onecomp(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);
ADVANFUNCS* pmx_advan_onecomp_depot(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);
ADVANFUNCS* pmx_advan_twocomp(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);
ADVANFUNCS* pmx_advan_threecomp(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);

/* eigensystem analytic advan methods */
ADVANFUNCS* pmx_advan_eigen(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);
ADVANFUNCS* pmx_advan_eigen_threecomp(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);
ADVANFUNCS* pmx_advan_eigen_twocomp(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);
ADVANFUNCS* pmx_advan_eigen_onecomp_absorb(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);

/* differential equation advan methods */
ADVANFUNCS* pmx_advan_diffeqn_libgsl(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);
ADVANFUNCS* pmx_advan_diffeqn_test(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);

/*---------------------------------------------------------------------*/
/* basic types */
/*---------------------------------------------------------------------*/
typedef struct {
	double objfn;
	enum {
		OBJFN_INVALID = 0,
		OBJFN_EVALUATE,
		OBJFN_CURRENT,
		OBJFN_FINAL,
	} type;
	int nparam;
	int neval;
} PMXRESULT;

/*---------------------------------------------------------------------*/
/* OPENPMX */
/*---------------------------------------------------------------------*/
typedef struct {
	/* options and config */
	const char* filename;
	int nthread;

	/* the problem structure to analyse */
	const DATACONFIG data;
	const ADVANCONFIG advan;

	/* the population paramaters and results of last run */
	struct {
		double lower, value, upper;
		enum {
			THETA_INVALID = 0,
			FIXED,
			ESTIMATE,
		} type;
	} theta[OPENPMX_THETA_MAX];
	struct {
		enum {
			OMEGA_INVALID = 0,
			OMEGA_DIAG,
			OMEGA_BLOCK,
			OMEGA_SAME,
		} type;
		int ndim;
		double values[OPENPMX_OMEGABLOCKSIZE_MAX];
	} omega[OPENPMX_OMEGABLOCK_MAX];
	double sigma[OPENPMX_SIGMA_MAX];

	/* output */
	PMXRESULT result;

	/* internal use */
	struct PMXSTATE* state;
} OPENPMX;

void pmx_cleanup(OPENPMX* openpmx);

OPENPMX pmx_copy(const OPENPMX* const openpmx);
void pmx_copy_popparam(OPENPMX* dest, const OPENPMX* const src);

/*---------------------------------------------------------------------*/
/* prediction */
/*---------------------------------------------------------------------*/
void pmx_predict(OPENPMX* openpmx);
void pmx_predict_pred(OPENPMX* openpmx);

/*---------------------------------------------------------------------*/
/* simulation */
/*---------------------------------------------------------------------*/
typedef struct {
	unsigned long seed;
} SIMCONFIG;

void pmx_resample(OPENPMX* pmx, const SIMCONFIG* const simconfig);
void pmx_simulate(OPENPMX* pmx, const SIMCONFIG* const simconfig);

/*---------------------------------------------------------------------*/
/* evaluation */
/*---------------------------------------------------------------------*/
typedef struct {
	double gradient_step;
	double step_initial;
	double step_refine;
	double step_final;
	int maxeval;

	bool icov_resample;
} STAGE1CONFIG;

void pmx_evaluate(OPENPMX* pmx, STAGE1CONFIG* const evalconfig);

/*---------------------------------------------------------------------*/
/* estimation */
/*---------------------------------------------------------------------*/
typedef struct {
	STAGE1CONFIG stage1;

	double step_initial;
	double step_refine;
	double step_final;
	int maxeval;
	double dobjfn;

	bool details;
	bool verbose;
} ESTIMCONFIG;

void pmx_estimate(OPENPMX* pmx, ESTIMCONFIG* const estimconfig);
void pmx_fastestimate(OPENPMX* pmx, ESTIMCONFIG* const estimconfig);

/*---------------------------------------------------------------------*/
/* tables */
/*---------------------------------------------------------------------*/
typedef struct {
	FILE* stream;
	const char* filename;
	const char* name;
	bool firstonly;
} TABLECONFIG;

void pmx_table(OPENPMX* pmx, const char* fields, const TABLECONFIG* const tableconfig);

/*---------------------------------------------------------------------*/
/* reload */
/*---------------------------------------------------------------------*/
typedef struct {
	const char* filename;
	const bool force;
	const bool optional;
	const bool preserve;
	const bool silent;
} RELOADCONFIG;

void pmx_reload_popparam(OPENPMX* dest, RELOADCONFIG* args);

/*---------------------------------------------------------------------*/
/* utility set functions */
/*---------------------------------------------------------------------*/
void pmx_set_theta(OPENPMX* dest, 
				   const int index,
				   typeof(((OPENPMX){0}).theta[0])* theta);
				   
/*---------------------------------------------------------------------*/
/* specialized model equations */
/*---------------------------------------------------------------------*/
typedef struct {
	/* input */
	const double age;
	const double weight;
	const double height;
	const bool male;
	
	/* output */
	double LBM;
	double V1, V2, V3, CL, Q2, Q3, KE0;
	double k10, k12, k21, k13, k31;
} SCHNIDER_PROPOFOL_CONFIG;

void pmx_model_schnider_propofol(SCHNIDER_PROPOFOL_CONFIG* schniderconfig);

#ifdef __cplusplus
}
#endif

/* OPENPMX_H */
#endif
