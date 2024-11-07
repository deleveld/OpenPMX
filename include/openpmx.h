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

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------------------------------------------------------*/
#define OPENPMX_VERSION_MAJOR			0
#define OPENPMX_VERSION_MINOR			0
#define OPENPMX_VERSION_RELEASE			1

#define OPENPMX_THETA_MAX				64
#define OPENPMX_OMEGABLOCKSIZE_MAX		64
#define OPENPMX_OMEGABLOCK_MAX			5
#define OPENPMX_OMEGA_MAX				16
#define OPENPMX_SIGMA_MAX				8

#define OPENPMX_FIELDS_MAX				64
#define OPENPMX_FIELDNAME_MAX			64

#define OPENPMX_STATE_MAX				32
#define OPENPMX_IMODEL_MAX				64
#define OPENPMX_PREDICTVARS_MAX			64
#define OPENPMX_SIMULINFUSION_MAX		16

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
typedef struct {
	/* info about the record being processed */
	const int initcount;
	const double statetime;
	const double recordtime;
	const RECORD* const record;

	/* may be initialized in IMODEL_INIT */
	double* const state;
	double* const amtlag;
	double* const bioavail;
	double _inittime;
} ADVANSTATE;

void pmx_advan_inittime(ADVANSTATE* const advanstate, const double t);

typedef	struct {
	const double* const theta;
	const int ntheta;
	const double* const eta;
	const int nomega;
	const double* const sigma;
	const int nsigma;
	const int nstate;
} POPPARAM;

/* must be initialized in IMODEL_INIT */
typedef struct IMODEL IMODEL;
typedef void (*IMODEL_INIT)(IMODEL* const imodel,
							ADVANSTATE* const advanstate,
							const POPPARAM* const popparam);

typedef struct {
	const double* const state;
	const RECORD* const record;
} PREDICTSTATE;

typedef struct PREDICTVARS PREDICTVARS;
typedef	double (*IMODEL_PREDICT)(const IMODEL* const imodel,
								 const PREDICTSTATE* const predictstate,
								 const POPPARAM* const popparam,
								 const double* const err,
								 PREDICTVARS* predparams);

typedef void (*ADVAN_DIFFEQN)(double DADT[],
							  const IMODEL* const imodel,
							  const RECORD* const record,
							  const double* state,
							  const POPPARAM* const popparam,
							  const double T);

/* unified arguments for all advancer functions */
typedef struct ADVANCONFIG ADVANCONFIG;
typedef struct ADVANFUNCS ADVANFUNCS;
typedef	ADVANFUNCS* (*ADVANMETHOD)(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);

typedef struct ADVANCONFIG {
	const IMODEL_INIT init;
	const IMODEL_PREDICT predict;
	const int firstonly;
	const int predictall;

	const STRUCTINFO imodelfields;
	const STRUCTINFO predictfields;
	const ADVANMETHOD method;
	const int nstate;
	const ADVAN_DIFFEQN diffeqn;

	/* extra arguments */
	union {
		/* for diffeqn solvers */
		struct {
			const double hstart;
			const double abstol;
			const double reltol;
			const char* steptype;
		} diffeqn;
		/* for wrapper solvers */
		struct {
			const ADVANMETHOD method;
		} inner;
	} args;
} ADVANCONFIG;

/* advancer methods */
ADVANFUNCS* pmx_advan_pred(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);
ADVANFUNCS* pmx_advan_onecomp(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);
ADVANFUNCS* pmx_advan_onecomp_depot(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);
ADVANFUNCS* pmx_advan_twocomp(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);
ADVANFUNCS* pmx_advan_threecomp(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);

ADVANFUNCS* pmx_advan_diffeqn_libgsl(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);

ADVANFUNCS* pmx_advan_diffeqn_test(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);
ADVANFUNCS* pmx_advan_wrapper(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig);

/*---------------------------------------------------------------------*/
/* basic types */
/*---------------------------------------------------------------------*/
typedef struct {
	double lower, value, upper;
	enum {
		THETA_INVALID = 0,
		FIXED = 1,
		ESTIMATE = 2,
	} type;
} THETA;

typedef struct {
	enum {
		OMEGA_INVALID = 0,
		OMEGA_DIAG,
		OMEGA_BLOCK,
		OMEGA_SAME,
	} type;
	int ndim;
	double values[OPENPMX_OMEGABLOCKSIZE_MAX];
} OMEGA;

typedef struct {
	double objfn;
	enum {
		OBJFN_INVALID = 0,
		OBJFN_PREDICT,
		OBJFN_EVALUATE,
		OBJFN_INITIAL,
		OBJFN_CURRENT,
		OBJFN_FINAL,
	} type;
	int nfunc;
} PMXRESULT;

/*---------------------------------------------------------------------*/
/* OPENPMX */
/*---------------------------------------------------------------------*/
typedef struct {
	/* options and config */
	const char* filename;
	int nthread;
	bool _offset1;
	bool verbose;
	bool brief;

	/* the problem structure to analyse */
	const DATACONFIG data;
	const ADVANCONFIG advan;

	/* the population paramaters and results of last run */
	THETA theta[OPENPMX_THETA_MAX];
	OMEGA omega[OPENPMX_OMEGABLOCK_MAX];
	double sigma[OPENPMX_SIGMA_MAX];
	PMXRESULT result;

	/* internal use */
	struct PMXSTATE* state;
} OPENPMX;

void pmx_cleanup(OPENPMX* openpmx);

void pmx_predict(OPENPMX* openpmx);
void pmx_predict_pred(OPENPMX* openpmx);

OPENPMX pmx_copy(const OPENPMX* const openpmx);
void pmx_copy_popparams(OPENPMX* dest, const OPENPMX* const src);

/*---------------------------------------------------------------------*/
/* estimation and simulation */
/*---------------------------------------------------------------------*/
/* for stage 1 (posthoc) estimation */
typedef struct {
	double gradient_step;
	double step_initial;
	double step_refine;
	double step_final;
	bool omit_icov_resample;
	int maxeval;
} STAGE1CONFIG;

/* for OBJFN minimization */
typedef struct {
	double step_initial;
	double step_refine;
	double step_final;
	int maxeval;
	bool fast;
} OPTIMCONFIG;

/* for simulations */
typedef struct {
	unsigned long seed;
} SIMCONFIG;

/* for estimation */
typedef struct {
	STAGE1CONFIG stage1;
	OPTIMCONFIG optim;
} ESTIMCONFIG;

void pmx_evaluate(OPENPMX* pmx, STAGE1CONFIG* const evalconfig);
void pmx_estimate(OPENPMX* pmx, ESTIMCONFIG* const estimconfig);
void pmx_fastestimate(OPENPMX* pmx, ESTIMCONFIG* const estimconfig);

void pmx_resample(OPENPMX* pmx, const SIMCONFIG* const simconfig);
void pmx_simulate(OPENPMX* pmx, const SIMCONFIG* const simconfig);

/*---------------------------------------------------------------------*/
/* tables */
/*---------------------------------------------------------------------*/
typedef struct {
	const char* filename;
	const char* name;
	bool firstonly;
} TABLECONFIG;

void pmx_table(OPENPMX* pmx, const char* fields, const TABLECONFIG* const tableconfig);

#ifdef __cplusplus
}
#endif

/*---------------------------------------------------------------------*/
/* C++ interface */
/*---------------------------------------------------------------------*/
#ifdef __cplusplus
namespace openpmx {
	class PPMX
	{
		private:
			::OPENPMX* _pmx;

			PPMX() = delete;
			PPMX(const PPMX&) = delete;

		public:
			PPMX(::OPENPMX& thepmx)
				:_pmx(&thepmx)
			{
				assert(_pmx->state == 0);
			};
			~PPMX()
			{
				::pmx_cleanup(_pmx);
			};

			void estimate(void)
			{
				::pmx_estimate(_pmx, 0);
			};
			void estimate(ESTIMCONFIG& estimconfig)
			{
				::pmx_estimate(_pmx, &estimconfig);
			};
			void evaluate(void)
			{
				::pmx_evaluate(_pmx, 0);
			};
			void evaluate(STAGE1CONFIG& stage1config)
			{
				::pmx_evaluate(_pmx, &stage1config);
			};
			void simulate(void)
			{
				::pmx_simulate(_pmx, 0);
			};
			void simulate(const SIMCONFIG& simconfig)
			{
				::pmx_simulate(_pmx, &simconfig);
			};
			void table(const char* fields)
			{
				::pmx_table(_pmx, fields, 0);
			};
			void table(const char* fields, const TABLECONFIG& tableconfig)
			{
				::pmx_table(_pmx, fields, &tableconfig);
			};

	};
};
#endif

/*---------------------------------------------------------------------*/
/* models */
/*---------------------------------------------------------------------*/
int remifentanil_eleveld(const double AGE,
						 const double WGT, 
						 const double HGT, 
						 const double M1F2,
						 double* const _V1, 
						 double* const _V2, 
						 double* const _V3, 
						 double* const _CL, 
						 double* const _Q2, 
						 double* const _Q3,
						 double* const _KE0);

/* OPENPMX_H */
#endif
