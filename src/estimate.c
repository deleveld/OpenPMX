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
 
/// This file does the outer (stage 2) estimation.
///
/// OpenPMX estimates population models in two nested stages. Stage 1
/// (see stage1.c) finds, for fixed population parameters, the individual
/// random effects (etas) that minimize each individual's objective
/// function contribution, and evaluates the first-order conditional
/// estimation with interaction (FOCEI) approximation of the likelihood.
/// Stage 2, implemented here, searches over the population parameters
/// (theta, omega and sigma) to minimize the total objective function
/// summed over all individuals. Every objective function evaluation of
/// stage 2 therefore runs stage 1 for all individuals, in parallel over
/// the available threads.
///
/// The public entry points are pmx_estimate() and pmx_evaluate().
/// Evaluation is estimation with a single objective function evaluation:
/// the model is not changed, but the individual etas, the objective
/// function, the predictions and the output files are produced exactly
/// as they would be for an estimation.
///
/// Parameter encoding
///
/// The population parameters are not optimized directly. They are encoded
/// (see encode.h) as an unconstrained vector of offsets relative to the
/// current model, so an all-zero vector reproduces the current model.
/// This puts parameters of very different scales on a common footing.
/// The step sizes used by the optimizer are therefore in encoded units,
/// not in the units of the parameters. Only the parameters that are
/// actually estimated are part of the vector. If there are none, no
/// optimization is done after the initial evaluation.
///
/// Optimization
///
/// The objective function is minimized with BOBYQA, a derivative-free
/// bounded quadratic-model trust region optimizer. Estimation runs in
/// phases, controlled by the estimate options (maxeval, step_initial,
/// step_refine, nsig, dobjfn):
///
///  1. Initial evaluation. The model as given is evaluated once, which
///     sets the starting objective function and the starting etas. A
///     warning is given if it is not finite.
///  2. Initial optimization. BOBYQA runs with a trust region radius
///     from step_initial down to step_refine. This makes large changes
///     to the parameters and does the bulk of the work.
///  3. Final step size. The smallest trust region radius that can still
///     resolve nsig significant digits in the parameters is found by
///     halving a radius, starting at half of step_refine, until
///     perturbing each encoded parameter by plus or minus that radius
///     leaves the decoded parameters agreeing to more than nsig digits.
///     The search stops at 1e-8.
///  4. Refinement. BOBYQA is restarted from the best model so far with a
///     radius from step_refine down to the final step size. After each
///     restart the number of significant digits (nsig) shared by the
///     parameters before and after the restart is reported, as the
///     minimum over all theta, sigma and the lower triangle of omega.
///     Restarts continue until nsig exceeds the requested nsig
///     (converged), until the change in the objective function between
///     restarts is smaller than dobjfn (stalled), until BOBYQA reports
///     an error, or until the evaluations allowed by maxeval are used up.
///
/// The result is marked as final only if the refinement converged,
/// otherwise it is the best model found so far. BOBYQA errors are
/// reported as warnings and do not discard the progress made.
///
/// Best model so far
///
/// The best model is tracked separately from the model under test and is
/// returned to the caller, so the result is always the best model seen
/// and not the last model evaluated. Only a finite objective function
/// is ever accepted as best: NaN and infinite values are rejected, and
/// a non-finite starting objective function is replaced by the first
/// finite one. When the best model improves, the etas of every individual
/// are saved and used as starting values for the next evaluations.
///
/// Output
///
/// Progress is written to the screen and, if the model has a filename, to
/// the out file. By default a single line is updated in place on a
/// terminal, the details and verbose options print every improvement or
/// every evaluation instead. At the start a header with the version,
/// the parallel configuration and the data is written, followed by a
/// data checkout (see checkout.c and ievaluate.c) that detects errors
/// before any time is spent on estimation. During estimation each
/// improvement is appended to the ext file. At the end the phi file
/// (individual parameters and etas) and the yhat file (pred, yhat,
/// yhatvar and the state) are written and the ext file gets its trailer.

#include "openpmx.h"
#include "githash.h"
#include "omegainfo.h"
#include "ievaluate.h"
#include "encode.h"
#include "stage1.h"
#include "defines.h"
#include "checkout.h"
#include "print.h"
#include "options.h"
#include "pmxstate.h"
#include "predict.h"
#include "advan/advan.h"
#include "bobyqa/bobyqa.h"
#include "utils/c22.h"
#include "utils/various.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <unistd.h>

#include "buildflags.h"

/*--------------------------------------------------------------------*/

typedef struct {
	IDATA* const idata;
	const ADVANFUNCS* const advanfuncs;
	ENCODE test;
	const OPTIONS* const options;
	struct {
		POPMODEL* model;
		IDATAETAS etas;
	} best;
	struct timespec begin;
	FILE* outstream;
	FILE* extstream;
} STAGE2_PARAMS;

static double get_timestamp(const STAGE2_PARAMS* const params)
{
	struct timespec now;
	clock_gettime(CLOCK_MONOTONIC, &now);
	return timespec_time_difference(&params->begin, &now) / 1000.;
}

/*--------------------------------------------------------------------*/
/* progress reporting: when not detailed or verbose, the current line is
 * overwritten in place on the terminal */

static bool quiet(const OPTIONS* const options)
{
	return !options->estimate.details && !options->estimate.verbose;
}

static void progress_open(const OPTIONS* const options)
{
	if (quiet(options))
		printf("\n");
}

static void progress_rewind(const OPTIONS* const options)
{
	if (quiet(options) && isatty(fileno(stdout)))
		printf("\033[F");
}

static void print_model_evaluation(const STAGE2_PARAMS* const params)
{
	let popmodel = &params->test.popmodel;
	let options = params->options;
	popmodel_eval_information(popmodel,
							  get_timestamp(params),
							  idata_ineval(params->idata, false),
							  !quiet(options),
							  params->outstream,
							  params->extstream,
							  0);
}

static void update_best_imodel(STAGE2_PARAMS* const params)
{
	let popmodel = &params->test.popmodel;
	let options = params->options;
	var best = params->best.model;

	/* update the best estimation and its objective function and function evaluations so far */
	/* only finite objective functions are ever accepted as best. A 
	 * non-finite best (only possible from the very first evaluation) is
	 * replaced by any finite objective function. */
	let objfn = popmodel->result.objfn;
	bool print_model = false;
	if (isfinite(objfn) &&
		(!isfinite(best->result.objfn) || objfn < best->result.objfn)) {
/// Each time the objective function improves during estimation the eta
/// values are saved. They will be used as initial values for subsequent
/// optimization.
		*best = *popmodel;
		idata_etas_copy(&params->best.etas, params->idata);
		print_model = true;
	}

	/* even if we didnt improve, update the number of evaluations of the best model */
	best->result.neval = popmodel->result.neval;

	/* update the user */
	if (options->estimate.verbose)
		print_model = true;
	if (print_model) {
		progress_rewind(options);
		print_model_evaluation(params);
	}
}

static void encode_evaluate(ENCODE* const test,
							IDATA* const idata,
							const ADVANFUNCS* const advanfuncs,
							const OPTIONS* const options,
							const bool first_eval)
{
	var popmodel = &test->popmodel;
	let omegainfo = &test->omegainfo;
	let nonzero = &omegainfo->nonzero;
	var scatteroptions = (SCATTEROPTIONS) {
		.stage1_order = true,
		.stage1_evaluate = first_eval,
	};
	scatter_threads(idata, advanfuncs, popmodel, nonzero, options, &scatteroptions, stage1_thread);

	popmodel->result.objfn = idata_objfn(idata, omegainfo->omega_nonzero_lndet);
	popmodel->result.type = OBJFN_CURRENT;
	popmodel->result.neval += 1;
	popmodel->result.nsig = 0.;
}

static double stage2_objfn(const long int _xlength,
						   const double _x[static _xlength],
						   void * const data)
{
	/* decode the vector in the right way for the advan and make the popmodel ready to test */
	let params = (STAGE2_PARAMS*) data;
	let idata = params->idata;
	assert(_xlength == params->test.nparam);
	encode_untransform(&params->test, _x);

	/* do the actual test, this sets the objfn */
	let advanfuncs = params->advanfuncs;
	let options = params->options;
	let popmodel = &params->test.popmodel;
	idata_etas_write(idata, &params->best.etas);
	encode_evaluate(&params->test, idata, advanfuncs, options, false);

	/* update best imodel and inform the user if we improve */
	update_best_imodel(params);

	return popmodel->result.objfn;
}

/*--------------------------------------------------------------------*/
/* BOBYQA */

static bool bobyqa_error(const int retcode, const char* phase, FILE* stream)
{
	if (retcode == BOBYQA_SUCCESS)
		return false;

	static const struct {
		int code;
		const char* msg;
	} errors[] = {
		{ BOBYQA_BAD_NPT,              "NPT is not in the required interval"             },
		{ BOBYQA_TOO_CLOSE,            "insufficient space between the bounds"           },
		{ BOBYQA_ROUNDING_ERRORS,      "too much cancellation in a denominator"          },
		{ BOBYQA_TOO_MANY_EVALUATIONS, "maximum number of function evaluations exceeded" },
		{ BOBYQA_STEP_FAILED,          "a trust region step has failed to reduce Q"      },
	};

	const char* errmsg = "unknown";
	forarray(i, errors) {
		if (errors[i].code == retcode) {
			errmsg = errors[i].msg;
			break;
		}
	}

	warning(stream, "%s BOBYQA error %i: %s\n", phase, retcode, errmsg);
	return true;
}

static long bobyqa_npt(const long n)
{
	return 2*n+1;				/* reccommended */
//	return n+2;					/* minimum */
//	return (n+1)*(n+2)/2;		/* maximum */
}

/* run one BOBYQA optimization starting from the current encoder
 * offsets (x = 0). Returns false if BOBYQA reported an error. */
static bool bobyqa_stage(STAGE2_PARAMS* const params,
						 const char* const phase,
						 const double rhobeg,
						 const double rhoend,
						 const int neval)
{
	let n = params->test.nparam;
	let npt = bobyqa_npt(n);
	let wsize = (npt+5)*(npt+n)+3*n*(n+5)/2 + 10; /* a little bit extra room to be sure */

	var x = mallocvar(double, n);
	var lower = mallocvar(double, n);
	var upper = mallocvar(double, n);
	var w = mallocvar(double, wsize);
	forcount(i, n) {
		x[i] = 0.;
		lower[i] = -1e6;	// -DBL_MAX;
		upper[i] = 1e6;		// DBL_MAX;
	}

	let iprint = 0;
	let retcode = bobyqa(n, npt,
						 stage2_objfn, (void*)params,
						 x, lower, upper,
						 rhobeg, rhoend,
						 iprint, neval, w);
	free(x);
	free(lower);
	free(upper);
	free(w);
	return !bobyqa_error(retcode, phase, params->outstream);
}

/*--------------------------------------------------------------------*/
/* number of significant digits */

static double compute_nsig_param(const double old_val, const double new_val)
{
	let delta = fabs(new_val - old_val);
	if (delta == 0.)
		return DBL_MAX;
	let scale = fmax(fabs(new_val), fabs(old_val));
	return fmax(0., -log10(delta / scale));
}

static double nsig_array(const double* const old_val,
						 const double* const new_val,
						 const long n)
{
	var min_nsig = DBL_MAX;
	forcount(i, n)
		min_nsig = fmin(min_nsig, compute_nsig_param(old_val[i], new_val[i]));
	return min_nsig;
}

static double iteration_nsig(const POPMODEL* const old_popmodel,
							 const POPMODEL* const popmodel)
{
	var min_nsig = fmin(nsig_array(old_popmodel->theta, popmodel->theta, popmodel->ntheta),
						nsig_array(old_popmodel->sigma, popmodel->sigma, popmodel->nsigma));

	/* lower triangle of omega */
	forcount(r, popmodel->nomega)
		min_nsig = fmin(min_nsig, nsig_array(old_popmodel->omega[r], popmodel->omega[r], r+1));

	return min_nsig;
}

static double converged_nsig(const POPMODEL* const popmodel, const double rhoend)
{
	var encoder = encode_init(popmodel);
	encode_transform(&encoder, popmodel);
	
	double x[encoder.nparam];
	memset(x, 0, sizeof(x));

	var nsig = DBL_MAX;
	forcount(i, encoder.nparam) {
		for (int sign = 1; sign >= -1; sign -= 2) {
			x[i] = sign * rhoend;
			encode_untransform(&encoder, x);
			nsig = fmin(nsig, iteration_nsig(popmodel, &encoder.popmodel));
		}
		x[i] = 0.;
	}
	return nsig;
}

/*--------------------------------------------------------------------*/
/* stage 2 optimization */

/// After initial optimization, a smaller search space is used from 
/// rho_refine to rho_final. This is repeated until the change in 
/// parameters between restarts is greater than nsig or objective 
/// function between restarts is less than dobjfn.
static bool refine(STAGE2_PARAMS* const params,
				   const double rhobeg,
				   const double rhoend,
				   int neval)
{
	let options = params->options;
	var best = params->best.model;
	var lastobjfn = best->result.objfn;

	for (;;) {
		if (neval <= bobyqa_npt(params->test.nparam)) {
			info(params->outstream, "optim too few iterations (%i) for refine\n", neval);
			return false;
		}

		/* set the encoder offsets and get ready to estimate again */
		let last_best = *best;
		encode_transform(&params->test, best);
		progress_open(options);

		/* setting etas to zero at each restart seems to reduce the 
		 * risk of local minimums */
		let idata = params->idata;
		idata_etas_reset(&params->best.etas, idata);

		let ok = bobyqa_stage(params, "stage2 refine", rhobeg, rhoend, neval);
		neval = options->estimate.maxeval - best->result.neval;
		if (!ok)
			return false;

		/* inform about the progress */
		char nsig_msg[16] = "(*)";
		let nsig = iteration_nsig(&last_best, best);
		if (nsig != DBL_MAX) {
			snprintf(nsig_msg, sizeof(nsig_msg), "%.1f", nsig);
			best->result.nsig = nsig;
		}
		let dobjfn = best->result.objfn - lastobjfn;
		progress_rewind(options);
		info(params->outstream,
			 "time %.3f neval %i objfn %.6f nsig %s dobjfn %.6f\n",
			 get_timestamp(params),
			 best->result.neval,
			 best->result.objfn,
			 nsig_msg,
			 dobjfn);

		/* check if we have converged to nsig */
		if (nsig > options->estimate.nsig)
			return true;

		/* check if we have no progress */
		lastobjfn = best->result.objfn;
		if (fabs(dobjfn) < fabs(options->estimate.dobjfn)) {
			info(params->outstream, "optim stalled nsig %s dobjfn %.6f\n", nsig_msg, dobjfn);
			return false;
		}
	}
}

static bool focei(STAGE2_PARAMS* const params)
{
	assert(params);
	let options = params->options;
	var best = params->best.model;

	if (params->test.nparam == 0) {
		info(params->outstream, "optim no parameters to estimate\n");
		return true;
	}

	let maxeval = options->estimate.maxeval;
	let step_initial = options->estimate.step_initial;
	let step_refine = options->estimate.step_refine;

/// Model estimation begins with an intial optimization with large 
/// changes to paramater values with BOBYQA rho values from 
/// step_initial, stopping at step_refine. 
	progress_open(options);
	/* warn on error, but try refine anyway */
	bobyqa_stage(params, "stage2 initial", step_initial, step_refine, maxeval);
	let neval = maxeval - best->result.neval;

	/* after initialization run find the final step size consistent with
	 * getting at least nsig */
	var step_final = step_refine / 2.;
	while (step_final >= 1e-8 &&
		   converged_nsig(best, step_final) <= options->estimate.nsig)
		step_final /= 2.;
	info(params->outstream, "optim rho %g %g\n", step_refine, step_final);

	let converged = neval > 1 && refine(params, step_refine, step_final, neval);
	if (converged) {
		best->result.nsig = converged_nsig(best, step_final);
//		print_converged_nsig(best, options->estimate.step_final);
		info(params->outstream, "optim final nsig %.1f\n", best->result.nsig);
	}
	return converged;
}

static void test_initial_objfn(STAGE2_PARAMS* const params)
{
	let idata = params->idata;
	let advanfuncs = params->advanfuncs;
	let options = params->options;
	let popmodel = &params->test.popmodel;

	/* start from best eta so far, is this going to be 0? */
	idata_etas_write(idata, &params->best.etas);

	/* very first evaluation */
	encode_evaluate(&params->test, idata, advanfuncs, options, true);
	info(params->outstream,
		 "time %.3f neval %i objfn %f\n",
		 get_timestamp(params),
		 popmodel->result.neval,
		 popmodel->result.objfn);
	if (!isfinite(popmodel->result.objfn))
		warning(params->outstream, "objective function not finite\n");
	
	/* save the popmodel, it shouldnt actually have been changed. It is 
	 * kept even if not finite, but will be replaced by the first finite 
	 * objective function. */
	*params->best.model = *popmodel;
	idata_etas_copy(&params->best.etas, idata);
}

static void focei_popmodel_stage2(STAGE2_PARAMS* const params)
{
	let idata = params->idata;
	let options = params->options;

	/* cleanup previous runs */
	/* at each estimate or evaluate, the pred and state gets reset to zero */
	let ndata = idata->ndata;
	var firstindivid = &idata->individ[0];
	memset(firstindivid->pred, 0, ndata * sizeof(double));
	memset(firstindivid->istate, 0, ndata * idata->nstate * sizeof(double));

	idata_free_simerr(idata);

	/* first run so we can set objective function and yhat. This also 
	 * saves the result in best. */
	test_initial_objfn(params);

	/* make sure 0 values reproduce the current model */
	var best = params->best.model;
	encode_transform(&params->test, best); /* this probably does nothing */

	/* call the underlying advan to optimize and find the best imodel */
	best->result.type = OBJFN_EVALUATE;
	if (options->estimate.maxeval > 1) {
		let fullconverge = focei(params); /* TODO: this should be renamed I think */
		if (fullconverge) 
			best->result.type = OBJFN_FINAL;

		/* at end we encode the best so far */
		encode_transform(&params->test, best);
	}
}

/*--------------------------------------------------------------------*/
/* setup, output and cleanup */

#if defined(OPENPMX_PARALLEL_PTHREADS)
#define PARALLEL_NAME "pthread"
#elif defined(OPENPMX_PARALLEL_OPENMP)
#define PARALLEL_NAME "openmp"
#elif defined(OPENPMX_PARALLEL_SINGLETHREAD)
#define PARALLEL_NAME "singlethread"
#else
#error no parallel processing defined
#endif

#if defined(OPENPMX_SERVER)
#define SERVER_NAME "server"
#elif defined(OPENPMX_NOTSERVER)
#define SERVER_NAME "notserver"
#else
#error no server option defined
#endif

static void outfile_header(FILE* f2,
						   const ADVANFUNCS* const advanfuncs,
						   const IDATA* const idata,
						   const OPTIONS* const options)
{
	assert(advanfuncs);
	assert(idata);
	assert(options);

	let dataconfig = advanfuncs->recordinfo.dataconfig;
	
	info(f2, "OpenPMX %i.%i.%i hash %s\n", OPENPMX_VERSION_MAJOR, OPENPMX_VERSION_MINOR, OPENPMX_VERSION_RELEASE, OPENPMX_GITHASH);
	info(f2, "config %s %i %s %s\n", 
		PARALLEL_NAME, options->nthread, SERVER_NAME, OPENPMX_INSTALL_PLATFORM);
	info(f2, "data records %i observations %i individuals %i\n", dataconfig->nrecords, idata->nobs, idata->nindivid);
	info(f2, "data table offset %s\n", dataconfig->_offset1 ? "true" : "false");

	advanfuncs->info(advanfuncs, stdout);
	if (f2)
		advanfuncs->info(advanfuncs, f2);
}

static FILE* open_or_die(FILE* const log, const char* const name, const char* const ext)
{
	char fname[PATH_MAX];
	snprintf(fname, sizeof(fname), "%s%s", name, ext);
	FILE* const f = fopen(fname, "w");
	if (!f)
		fatal(log, "%s: could not open file \"%s\"\n", __func__, fname);
	return f;
}

static STAGE2_PARAMS stage2_params_init(const char* filename,
										IDATA* const idata,
										const ADVANFUNCS* const advanfuncs,
										POPMODEL* const popmodel,
										const OPTIONS* const options)
{
	if (idata->nindivid <= 0)
		fatal(0, "%s: optim cannot estimate, no individuals\n", __func__);

	/* setup output logging */
	FILE* outstream = 0;
	FILE* extstream = 0;
	if (filename) {
		outstream = open_or_die(0, filename, OPENPMX_OUTFILE);
		extstream = open_or_die(outstream, filename, OPENPMX_EXTFILE);
	}

	struct timespec begin;
	clock_gettime(CLOCK_MONOTONIC, &begin);

	var params = (STAGE2_PARAMS) {
		.idata = idata,
		.advanfuncs = advanfuncs,
		.test = encode_init(popmodel),	/* makes a copy of popmodel to work with */
		.options = options,
		.best = {
			.model = popmodel,			/* best so far will be saved at the caller */
			.etas = idata_etas_alloc(idata),
		},
		.begin = begin,
		.outstream = outstream,
		.extstream = extstream,
	};

	/* save the current, but with objfn to invalid and nparam as well */
	*params.best.model = params.test.popmodel;

	/* make sure 0 values reproduce the current model */
	encode_transform(&params.test, params.best.model);

	return params;
}

static void stage2_params_cleanup(STAGE2_PARAMS* const params)
{
	if (params->outstream)
		fclose(params->outstream);
	if (params->extstream)
		fclose(params->extstream);
	idata_etas_free(&params->best.etas);
}

static void estimate_popmodel(const char* filename,
							  IDATA* const idata,
							  const ADVANFUNCS* const advanfuncs,
							  POPMODEL* const popmodel,
							  const OPTIONS* const options)
{
	let estimating = options->estimate.maxeval > 1;
	popmodel->result.type = OBJFN_CURRENT;

	/* setup the minimizer */
	var params = stage2_params_init(filename, idata, advanfuncs, popmodel, options);

/// At start of estimation the header of the ext file is written.
	/* do some logging */
	/* start extfile and other headers, rest will be saved during iterations */
	outfile_header(params.outstream, advanfuncs, idata, options);
	let _offset1 = advanfuncs->recordinfo.dataconfig->_offset1;
	if (params.extstream) 
		extfile_header(params.extstream, params.best.model, _offset1, true);
	if (estimating)
		info(params.outstream, "optim FOCEI BOBYQA\n");

	/* update results to screen and log */
	popmodel_information(params.outstream, popmodel, -DBL_MAX);

/// Before estimation a data checkout is done (see checkout.c and 
/// ievaluate.c) to detect various errors. 
	idata_ineval(idata, true);
	idata_checkout(idata, advanfuncs, popmodel, options, params.outstream);

	/* actually do the stage 2 estimation */
	/* best.model in params saves the best at the caller so we dont
	 * need to save it ourselves */
	focei_popmodel_stage2(&params);

	/* update results to screen and log */
	let timestamp = get_timestamp(&params);
	popmodel_information(params.outstream,
						 popmodel,
						 timestamp);

/// At the end of estimation the phi file is written and a trailer is 
/// put onto the ext file. Also the yhat file is written with prediction
/// varaiables including pred, yhat, yhatvar, and the state.
	/* update results in tables */
	if (filename) {
		table_phi_idata(filename, idata, _offset1);

		/* TODO: could this be optional, it might be burdensome for 
		 * large models because it forces another advancing throught the 
		 * records */
		idata_predict_pred(idata, advanfuncs, popmodel, options);
		table_yhat_idata(filename, idata, _offset1);

		if (params.extstream) {
			let ineval = idata_ineval(idata, false);
			extfile_trailer(params.extstream, popmodel, timestamp, ineval);
		}
	}
	
	/* cleanup */
	stage2_params_cleanup(&params);
}

static void run_popmodel(OPENPMX* const pmx, const OPTIONS* const options)
{
	var pstate = pmx->state;

	ERRCTX errctx = { 0 };
	var popmodel = popmodel_init(pmx->theta, pmx->omega, pmx->sigma, &errctx);
	if (errctx.len)
		fatal(0, "%s: %s", __func__, errctx.errmsg);

	estimate_popmodel(pmx->filename,
					  &pstate->idata,
					  pstate->advanfuncs,
					  &popmodel,
					  options);

	pmx_copy_popmodel(pmx, &popmodel);
}

/*--------------------------------------------------------------------*/
/* public interface */

void pmx_estimate(OPENPMX* pmx, ESTIMCONFIG* const estimate)
{
	pmx_ensure_state(pmx);
	var pstate = pmx->state;

	var options = options_init(pmx);
	if (estimate) {
		options.estimate = estimconfig_default(estimate);
		*estimate = options.estimate;
	}

	/* estimate destroys a previous covariance matrix */
	pmx_covariance_free(&pstate->covariance);
	pstate->covariance = (COVARIANCE){ 0 };

	run_popmodel(pmx, &options);
}

/// Evaluation is the same as estimation but with maxeval=1. Note that 
/// maxeval=0 would trigger the default number of evaluations.

void pmx_evaluate(OPENPMX* pmx, STAGE1CONFIG* const stage1)
{
	pmx_ensure_state(pmx);

	var options = options_init(pmx);
	if (stage1) {
		options.estimate.stage1 = stage1config_default(stage1);
		*stage1 = options.estimate.stage1;
	}
	options.estimate.maxeval = 1;

	run_popmodel(pmx, &options);
}
