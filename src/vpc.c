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

/*
 * Prediction-corrected Visual Predictive Check (pcVPC).
 *
 * Overview of pmx_vpc():
 *   1. Validate/default the config and set up a private copy of the model.
 *   2. Compute PRED (population prediction) for every record.
 *   3. Select the usable observations (DVTY filter, finite values, ...).
 *   4. Split observations into strata and equal-count bins along x.
 *   5. Compute observed (optionally prediction-corrected) percentiles per bin.
 *   6. Simulate nsim replicates, computing the same percentiles per bin.
 *   7. Write the summary and points files via temp files + rename, so a
 *      failure never leaves a half-written output behind.
 */

#include "openpmx.h"
#include "idata.h"
#include "popmodel.h"
#include "print.h"
#include "options.h"
#include "scatter.h"
#include "simulate.h"
#include "pmxstate.h"
#include "defines.h"
#include "advan/advan.h"
#include "utils/c22.h"
#include "utils/errctx.h"
#include "utils/filecontent.h"
#include "utils/various.h"

#include <ctype.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include <gsl/gsl_statistics.h>

/* rename() fails on Windows if the target exists, so replace explicitly.
 * Returns 0 on success like rename(). */
#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#define VPC_RENAME(from, to) \
	(MoveFileExA((from), (to), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH) ? 0 : -1)
#else
#define VPC_RENAME(from, to) rename((from), (to))
#endif

/* forward declaration (defined at the bottom of this file) */
void pmx_vpc_clean(const OPENPMX* const openpmx, const char* name);

/* x value paired with the index o of the observation it belongs to, so a
 * stratum can be sorted by x while remembering the original observation */
typedef struct { double x; int o; } VPC_XO;

/* qsort comparator for doubles, ascending. (x > y) - (x < y) avoids the
 * overflow/NaN pitfalls of returning a subtraction */
static int vpc_cmp_double(const void* a, const void* b)
{
	let x = *(const double*)a;
	let y = *(const double*)b;
	return (x > y) - (x < y);
}

/* qsort comparator for VPC_XO, ascending by x */
static int vpc_cmp_xo(const void* a, const void* b)
{
	let x = ((const VPC_XO*)a)->x;
	let y = ((const VPC_XO*)b)->x;
	return (x > y) - (x < y);
}

/* offset of a named field in the record, -1 if the field is absent.
 * Name matching is case insensitive. */
static int vpc_field_find(const STRUCTINFO* info, const char* name)
{
	forcount(i, OPENPMX_FIELDS_MAX)
		if (info->field[i].name[0] && strcasecmp(info->field[i].name, name) == 0)
			return info->field[i].offset;
	return -1;
}

/* like vpc_field_find() but a missing field is a fatal error */
static int vpc_field_offset(const STRUCTINFO* info, const char* name)
{
	let off = vpc_field_find(info, name);
	if (off < 0)
		fatal(0, "VPC field '%s' not found in record\n", name);
	return off;
}

/* copy the rows of filename whose first column != name into fout. The
 * header line is kept. *any is set if something was written. Returns
 * false if filename could not be read. */
static bool vpc_filter_to(FILE* fout, const char* filename, const char* name, bool* const any)
{
	/* read the whole file into a mutable buffer; failure (e.g. the file
	 * does not exist yet) is reported to the caller, not fatal */
	ERRCTX errctx = { 0 };
	var text = filecontent_alloc(filename, &errctx);
	if (errctx.len)
		return false;

	const size_t n = strlen(name);
	int lineno = 0;
	var ptr = text.mutptr;
	/* walk the buffer line by line, terminating each line in place */
	while (*ptr != '\0') {
		var next_line = strchr(ptr, '\n');
		if (next_line)
			*next_line = '\0';

		/* skip leading whitespace, then compare the first column with name.
		 * The column must end at whitespace/end of line so "A" does not
		 * match "AB". Line 0 is the header and is never removed. */
		const char* p = ptr;
		while (isspace((unsigned char)*p))
			p++;
		let match = lineno > 0 && n > 0 &&
					strncmp(p, name, n) == 0 &&
					(p[n] == '\0' || isspace((unsigned char)p[n]));
		/* keep non-matching, non-blank lines */
		if (!match && *p != '\0') {
			fprintf(fout, "%s\n", ptr);
			*any = true;
		}
		lineno++;

		if (!next_line)
			break;
		ptr = next_line + 1;
	}
	filecontent_free(&text);
	return true;
}

/* open a temporary file holding the existing rows of filename without
 * those of name, ready for appending. *empty is set if it has no
 * content yet. */
static FILE* vpc_open_tmp(const char* filename,
						  const char* name,
						  char* const tmp,
						  const size_t tmpsize,
						  bool* const empty)
{
	/* temp file lives next to the original so the later rename stays on
	 * the same filesystem */
	snprintf(tmp, tmpsize, "%s.tmp", filename);
	var f = fopen(tmp, "w");
	if (!f)
		fatal(0, "cannot open VPC output '%s'\n", tmp);
	/* carry over old rows (minus this VPC name). If the original does not
	 * exist, vpc_filter_to() writes nothing and the file counts as empty,
	 * so the caller knows to write a header */
	bool any = false;
	vpc_filter_to(f, filename, name, &any);
	*empty = !any;
	return f;
}

/* close a file and report whether all writes (including the final flush)
 * succeeded */
static bool vpc_close(FILE* f)
{
	bool ok = !ferror(f);
	if (fclose(f) != 0)
		ok = false;
	return ok;
}

/* write a space-separated number, or NA for NaN/Inf so the output stays
 * readable by R and similar tools */
static void vpc_num(FILE* out, const double v)
{
	if (isfinite(v))
		fprintf(out, " %.10g", v);
	else
		fprintf(out, " NA");
}

/* Bergstrand 2011 pcVPC: lb + (y-lb) * (PREDbin-lb) / (PRED-lb)
 * Rescales a value by the ratio of the bin's typical prediction to the
 * observation's own prediction, removing variability due to differing
 * covariates/doses. With no correction the value is returned unchanged. */
static inline double vpc_correct(const VPCCONFIG* const c,
								 const double y,
								 const double pred,
								 const double predbin)
{
	if (c->correction == VPC_NOCORRECTION)
		return y;
	let lb = c->lbound;
	return lb + (y - lb) * (predbin - lb) / (pred - lb);
}

/* percentiles (R type 7) of finite values within each bin.
 *   val   : one value per observation
 *   order : observation indices grouped by bin
 *   start : bin g owns order[start[g] .. start[g+1]-1] (nbin+1 entries)
 *   prob  : the nprob probabilities requested
 *   buf   : scratch space, at least as large as the biggest bin
 *   out   : nbin*nprob results, row g holds the percentiles of bin g.
 * A bin without finite values gets NaN for every percentile. */
static void vpc_percentiles(const double* const val,
							const int* const order,
							const int* const start,
							const int nbin,
							const double* const prob,
							const int nprob,
							double* const buf,
							double* const out)
{
	forcount(g, nbin) {
		/* gather the finite values of this bin into the scratch buffer */
		int m = 0;
		for (int j = start[g]; j < start[g + 1]; j++) {
			let v = val[order[j]];
			if (isfinite(v))
				buf[m++] = v;
		}
		if (m == 0) {
			forcount(k, nprob)
				out[g * nprob + k] = NAN;
			continue;
		}
		/* sort once, then read every requested quantile from it */
		qsort(buf, m, sizeof(double), vpc_cmp_double);
		forcount(k, nprob)
			out[g * nprob + k] = gsl_stats_quantile_from_sorted_data(buf, 1, m, prob[k]);
	}
}

/* run one simulation pass over all individuals using the ETAs and
 * residual errors currently stored in idata. Stage 1 ordering is off since
 * no estimation is done. Results end up in idata->individ[0].yhat / pred. */
static void vpc_simulate_pass(IDATA* const idata,
							  const ADVANFUNCS* const advanfuncs,
							  const POPMODEL* const popmodel,
							  const OPTIONS* const options)
{
	SCATTEROPTIONS scatteroptions = { 
		.stage1_order = false
	};
	scatter_threads(idata, advanfuncs, popmodel, 0, options, &scatteroptions, simulate_with_error_thread);
}

/// Prediction-corrected VPC. The dataset DV is not modified. Like
/// pmx_simulate(), the individual ETAs in the state are overwritten.
/// The results are appended to the VPC output and points files. The
/// first column is the VPC name and any previous VPC with the same name
/// is replaced. The new files are written to temporary files and renamed
/// over the old ones only after the VPC succeeded and was written fully.
void pmx_vpc(OPENPMX* source, VPCCONFIG* const cfg)
{
	info(0, "VPC begin\n");

	if (!cfg)
		fatal(0, "VPC config required\n");
	if (!source->filename || !*source->filename)
		fatal(0, "VPC needs a filename\n");

	/* ---- apply defaults / validate config ----
	 * nsim:  replicates (1000),
	 * nbins: bins per stratum (10),
	 * ci:    width of the simulated confidence interval (0.95),
	 * prob:  percentiles to report (5%, 50%, 95%) */
	let nsim = cfg->nsim > 0 ? cfg->nsim : 1000;
	let nbins = cfg->nbins > 0 ? cfg->nbins : 10;
	let ci = (cfg->ci > 0. && cfg->ci < 1.) ? cfg->ci : 0.95;
	double prob[OPENPMX_VPC_PROB_MAX] = { 0.05, 0.5, 0.95 };
	int nprob = 3;
	if (cfg->nprob > 0) {
		if (cfg->nprob > OPENPMX_VPC_PROB_MAX)
			fatal(0, "VPC nprob exceeds OPENPMX_VPC_PROB_MAX\n");
		nprob = cfg->nprob;
		memcpyvar(prob, cfg->prob, nprob);
		forcount(k, nprob)
			if (!(prob[k] >= 0. && prob[k] <= 1.))
				fatal(0, "VPC prob must be in [0,1]\n");
	}

	/* ---- private model copy and state ----
	 * Work on a copy so the caller's state is not touched (its ETAs are
	 * not overwritten), then build the population model, options and RNG */
	var pmx = pmx_copy(source);
	pmx.filename = source->filename;
	pmx_ensure_state(&pmx);
	var pstate = pmx.state;

	ERRCTX errctx = { 0 };
	var popmodel = popmodel_init(pmx.theta, pmx.omega, pmx.sigma, &errctx);
	if (errctx.len)
		fatal(0, "%s", errctx.errmsg);

	var options = options_init(&pmx);
	options.simulate = simconfig_default(&cfg->simulate);
	pmx_ensure_state_rng(&pmx, &options);

	/* ---- shortcuts into the data/state ---- */
	var idata = &pstate->idata;
	let advanfuncs = pstate->advanfuncs;
	let recordinfo = &advanfuncs->recordinfo;
	let dataconfig = recordinfo->dataconfig;
	let ndata = recordinfo->ndata;
	let nomega = idata->nomega;
	let nsigma = popmodel.nsigma;
	let recordsize = dataconfig->recordfields.size;
	/* record offsets of the x variable (default TIME) and the optional
	 * stratification variable (-1 = no stratification) */
	let offx = vpc_field_offset(&dataconfig->recordfields, cfg->xfield ? cfg->xfield : "TIME");
	let offs = cfg->stratfield ? vpc_field_offset(&dataconfig->recordfields, cfg->stratfield) : -1;

	/* DVTY filter field: the configured one, else DVTY or DVID if the
	 * dataset has such a field. The value is cfg->dvty.value. */
	var dvtyfield = cfg->dvty.field;
	if (!dvtyfield) {
		let fields = &dataconfig->recordfields;
		if (vpc_field_find(fields, "DVTY") >= 0)
			dvtyfield = "DVTY";
		else if (vpc_field_find(fields, "DVID") >= 0)
			dvtyfield = "DVID";
	}
	let offd = dvtyfield ? vpc_field_offset(&dataconfig->recordfields, dvtyfield) : -1;

	/* the name is the first column, it cannot be empty or have whitespace
	 * if name is missing use the dvty field and value */
	char autoname[128];
	var vpcname = (cfg->name && *cfg->name) ? cfg->name : "-";
	if (!(cfg->name && *cfg->name) && dvtyfield) {
		snprintf(autoname, sizeof(autoname), "%s%.10g", dvtyfield, cfg->dvty.value);
		vpcname = autoname;
	}
	/* whitespace would break the space-separated output format */
	for (var c = vpcname; *c; c++)
		if (isspace((unsigned char)*c))
			fatal(0, "VPC name '%s' cannot contain whitespace\n", vpcname);

	char dvtystr[32];	/* DVTY column value, NA if no DVTY filter */
	if (offd >= 0)
		snprintf(dvtystr, sizeof(dvtystr), "%.10g", cfg->dvty.value);
	else
		snprintf(dvtystr, sizeof(dvtystr), "NA");

	/* PRED: population prediction, ETA=0 and no residual error */
	forcount(i, idata->nindivid)
		memset(idata->individ[i].eta, 0, nomega * sizeof(double));
	memset(idata_alloc_simerr(idata), 0, ndata * nsigma * sizeof(double));
	vpc_simulate_pass(idata, advanfuncs, &popmodel, &options);
	var pred0 = mallocvar(double, ndata);
	memcpyvar(pred0, idata->individ[0].yhat, ndata);

	/* ---- observations used in the VPC ----
	 * Per-observation arrays (indexed 0..nobs-1):
	 *   oidx : record index in the dataset
	 *   ox   : x value        os   : stratum value
	 *   odv  : observed DV    opred: population prediction */
	var oidx = mallocvar(int, ndata);
	var ox = mallocvar(double, ndata);
	var os = mallocvar(double, ndata);
	var odv = mallocvar(double, ndata);
	var opred = mallocvar(double, ndata);
	int nobs = 0;
	forcount(i, ndata) {
		let ptr = (const char*)dataconfig->records + recordsize * i;
		/* only observation records (EVID == 0) */
		if (RECORDINFO_EVID(recordinfo, (const RECORD*)ptr) != 0)
			continue;
		/* only the selected DVTY */
		if (offd >= 0 && *(const double*)(ptr + offd) != cfg->dvty.value)
			continue;
		let dv = RECORDINFO_DV(recordinfo, (const RECORD*)ptr);
		let x = *(const double*)(ptr + offx);
		let p = pred0[i];
		let sv = (offs >= 0) ? *(const double*)(ptr + offs) : 0.;
		/* drop records with missing/non-finite values */
		if (!isfinite(dv) || !isfinite(x) || !isfinite(p) || !isfinite(sv))
			continue;
		/* the correction divides by (PRED - lbound), so PRED must be above it */
		if (cfg->correction != VPC_NOCORRECTION && p <= cfg->lbound)
			continue;
		oidx[nobs] = i;
		ox[nobs] = x;
		os[nobs] = sv;
		odv[nobs] = dv;
		opred[nobs] = p;
		nobs++;
	}
	if (nobs == 0)
		fatal(0, "VPC has no usable observations\n");

	/* strata: sort a copy of the stratum values and compact it in place
	 * to the unique values (ustrat[0..nstrat-1], ascending) */
	var ustrat = mallocvar(double, nobs);
	memcpyvar(ustrat, os, nobs);
	qsort(ustrat, nobs, sizeof(double), vpc_cmp_double);
	int nstrat = 0;
	forcount(i, nobs)
		if (nstrat == 0 || ustrat[i] != ustrat[nstrat - 1])
			ustrat[nstrat++] = ustrat[i];

	/* equal-count bins within each stratum, tied x never split, cuts snapped to nearest x change.
	 * Outputs: binof[o] = global bin index of observation o,
	 *          bstrat[b]/bnum[b] = stratum value / number within stratum of bin b,
	 *          bedlo[b]/bedhi[b] = bin edges (histogram limits) of bin b: the midpoint
	 *          between the last x of one bin and the first x of the next, so
	 *          neighbouring bins share an edge; the outer edges of a stratum are
	 *          its smallest and largest x */
	var binof = mallocvar(int, nobs);
	var bstrat = mallocvar(double, nobs);
	var bnum = mallocvar(int, nobs);
	var bedlo = mallocvar(double, nobs);
	var bedhi = mallocvar(double, nobs);
	var xo = mallocvar(VPC_XO, nobs);
	int nbin = 0;
	forcount(st, nstrat) {
		/* collect this stratum's observations and sort them by x */
		int n = 0;
		forcount(o, nobs)
			if (os[o] == ustrat[st])
				xo[n++] = (VPC_XO){ ox[o], o };
		qsort(xo, n, sizeof(VPC_XO), vpc_cmp_xo);
		/* never more bins than observations */
		let nb = n < nbins ? n : nbins;
		/* bin starts: ideal equal-count cut j*n/nb, moved to the nearest index
		 * where x changes so tied x are never split. Searching both ways (not
		 * only forward) keeps the highest-x observations from being absorbed
		 * into the previous bin when they are tied or sparse. */
		var cut = mallocvar(int, nb > 0 ? nb : 1);
		int ncut = 0;
		if (n > 0)
			cut[ncut++] = 0;
		for (int j = 1; j < nb; j++) {
			let ideal = (int)((long)j * n / nb);
			let lo = cut[ncut - 1] + 1;	/* a cut must lie after the previous one */
			int found = -1;
			/* expand outwards from the ideal index (backward kb, forward kf)
			 * until an index is found where x differs from its predecessor */
			for (int d = 0; found < 0 && (ideal - d >= lo || ideal + d < n); d++) {
				let kb = ideal - d;
				let kf = ideal + d;
				if (kb >= lo && kb < n && xo[kb].x != xo[kb - 1].x)
					found = kb;
				else if (kf >= lo && kf < n && xo[kf].x != xo[kf - 1].x)
					found = kf;
			}
			/* no valid cut left (too many ties): use fewer bins */
			if (found < 0)
				break;
			cut[ncut++] = found;
		}
		/* walk the sorted observations, opening a new bin at each cut and
		 * assigning every observation to the current bin */
		int local = 0, c = 0;
		forcount(k, n) {
			if (c < ncut && k == cut[c]) {
				/* shared edge with the previous bin of this stratum */
				let edge = (k == 0) ? xo[0].x : 0.5 * (xo[k - 1].x + xo[k].x);
				if (k > 0)
					bedhi[nbin - 1] = edge;
				bstrat[nbin] = ustrat[st];
				bnum[nbin] = ++local;
				bedlo[nbin] = edge;
				nbin++;
				c++;
			}
			binof[xo[k].o] = nbin - 1;
		}
		/* outer upper edge of the stratum */
		if (n > 0)
			bedhi[nbin - 1] = xo[n - 1].x;
		free(cut);
	}

	/* members of each bin: counting sort of observations by bin.
	 * Bin g owns order[start[g] .. start[g+1]-1] */
	var start = callocvar(int, nbin + 1);
	forcount(o, nobs)
		start[binof[o] + 1]++;			/* count per bin ... */
	forcount(g, nbin)
		start[g + 1] += start[g];		/* ... turned into start offsets */
	var fill = mallocvar(int, nbin);
	memcpyvar(fill, start, nbin);
	var order = mallocvar(int, nobs);
	forcount(o, nobs)
		order[fill[binof[o]]++] = o;

	/* bin summaries: median x, median PRED (for the correction).
	 * buf is shared scratch, sized for the largest user (a bin or nsim values) */
	var buf = mallocvar(double, (nobs > nsim ? nobs : nsim));
	var xmed = mallocvar(double, nbin);
	var predbin = mallocvar(double, nbin);
	forcount(g, nbin) {
		let n = start[g + 1] - start[g];
		forcount(j, n)
			buf[j] = ox[order[start[g] + j]];
		qsort(buf, n, sizeof(double), vpc_cmp_double);
		xmed[g] = gsl_stats_quantile_from_sorted_data(buf, 1, n, 0.5);
		forcount(j, n)
			buf[j] = opred[order[start[g] + j]];
		qsort(buf, n, sizeof(double), vpc_cmp_double);
		predbin[g] = gsl_stats_quantile_from_sorted_data(buf, 1, n, 0.5);
	}

	/* observed (corrected) values and percentiles. obsval is kept for the
	 * points file, val is reused as scratch for the simulations */
	var val = mallocvar(double, nobs);
	var obsval = mallocvar(double, nobs);
	forcount(o, nobs)
		obsval[o] = vpc_correct(cfg, odv[o], opred[o], predbin[binof[o]]);
	var obspct = mallocvar(double, (size_t)nbin * nprob);
	vpc_percentiles(obsval, order, start, nbin, prob, nprob, buf, obspct);

	struct timespec begin, now;
	clock_gettime(CLOCK_MONOTONIC, &begin);
	var stime = 0.;

	/* simulated percentiles for each replicate.
	 * simpct layout: [replicate][bin][prob] */
	var simpct = mallocvar(double, (size_t)nsim * nbin * nprob);
	forcount(s, nsim) {
		clock_gettime(CLOCK_MONOTONIC, &now);
		let t = timespec_time_difference(&begin, &now) / 1000.;
		if (t > stime + 1) {
			stime += 1.;
			info(0, "VPC sim %i\n", s + 1);
		} 
		/* draw new individual ETAs and residual errors, then simulate */
		idata_resample_eta(idata, &popmodel, pstate->rng);
		idata_resample_err(idata, &popmodel, pstate->rng);
		vpc_simulate_pass(idata, advanfuncs, &popmodel, &options);

		/* correct the simulated values exactly like the observed ones and
		 * take the per-bin percentiles of this replicate */
		let pred = idata->individ[0].pred;	/* observation with noise */
		forcount(o, nobs)
			val[o] = vpc_correct(cfg, pred[oidx[o]], opred[o], predbin[binof[o]]);
		vpc_percentiles(val, order, start, nbin, prob, nprob, buf, &simpct[(size_t)s * nbin * nprob]);
	}
	info(0, "VPC sim %i\n", nsim);

	/* everything succeeded. Build the new files next to the old ones
	 * (old rows of this name removed, new rows appended) and only
	 * rename them over the originals once both were written fully. */
	char fout_name[PATH_MAX], fpts_name[PATH_MAX];
	char tmpout[PATH_MAX + 8], tmppts[PATH_MAX + 8];
	bool out_empty, pts_empty;
	
	snprintf(fout_name, sizeof(fout_name), "%s" OPENPMX_VPC_OUT, pmx.filename);
	snprintf(fpts_name, sizeof(fpts_name), "%s" OPENPMX_VPC_POINTS, pmx.filename);
	var out = vpc_open_tmp(fout_name, vpcname, tmpout, sizeof(tmpout), &out_empty);
	var pts = vpc_open_tmp(fpts_name, vpcname, tmppts, sizeof(tmppts), &pts_empty);

	/* points file: one row per observation with its bin, raw DV, PRED and
	 * the (prediction-corrected) value. Header only if the file is new. */
	if (pts_empty)
		fprintf(pts, "NAME DVTY STRAT BIN X DV PRED PCDV\n");
	forcount(o, nobs) {
		fprintf(pts, "%s %s %.10g %i", vpcname, dvtystr, bstrat[binof[o]], bnum[binof[o]]);
		vpc_num(pts, ox[o]);
		vpc_num(pts, odv[o]);
		vpc_num(pts, opred[o]);
		vpc_num(pts, obsval[o]);
		fprintf(pts, "\n");
	}

	/* summary file: one row per bin with its bin edges (BINMIN/BINMAX), the observed
	 * percentiles and the confidence interval (lo/med/hi) of each
	 * percentile across the simulated replicates. Header only if new. */
	if (out_empty) {
		fprintf(out, "NAME DVTY STRAT BIN NOBS XMED BINMIN BINMAX");
		forcount(k, nprob) {
			char pname[32];
			snprintf(pname, sizeof(pname), "P%g", 100. * prob[k]);
			fprintf(out, " OBS_%s SIM_%s_LO SIM_%s_MED SIM_%s_HI", pname, pname, pname, pname);
		}
		fprintf(out, "\n");
	}
	/* tail probability on each side of the interval, e.g. 0.025 for ci=0.95 */
	let alpha = 0.5 * (1. - ci);
	forcount(g, nbin) {
		fprintf(out, "%s %s %.10g %i %i", vpcname, dvtystr, bstrat[g], bnum[g], start[g + 1] - start[g]);
		vpc_num(out, xmed[g]);
		vpc_num(out, bedlo[g]);
		vpc_num(out, bedhi[g]);
		forcount(k, nprob) {
			/* distribution of this bin's k-th percentile over the replicates */
			int m = 0;
			forcount(s, nsim) {
				let v = simpct[((size_t)s * nbin + g) * nprob + k];
				if (isfinite(v))
					buf[m++] = v;
			}
			double lo = NAN, med = NAN, hi = NAN;
			if (m > 0) {
				qsort(buf, m, sizeof(double), vpc_cmp_double);
				lo = gsl_stats_quantile_from_sorted_data(buf, 1, m, alpha);
				med = gsl_stats_quantile_from_sorted_data(buf, 1, m, 0.5);
				hi = gsl_stats_quantile_from_sorted_data(buf, 1, m, 1. - alpha);
			}
			vpc_num(out, obspct[g * nprob + k]);
			vpc_num(out, lo);
			vpc_num(out, med);
			vpc_num(out, hi);
		}
		fprintf(out, "\n");
	}

	/* commit: close both files, and only if both were written without
	 * error replace the originals. On any failure delete the temp files
	 * so the existing outputs stay intact. */
	let ok_out = vpc_close(out);
	let ok_pts = vpc_close(pts);
	if (!ok_out || !ok_pts) {
		remove(tmpout);
		remove(tmppts);
		fatal(0, "error writing VPC output\n");
	}
	if (VPC_RENAME(tmpout, fout_name) != 0 || VPC_RENAME(tmppts, fpts_name) != 0) {
		remove(tmpout);
		remove(tmppts);
		fatal(0, "cannot replace VPC output files\n");
	}

	/* ---- cleanup ---- */
	pmx_release_state(&pmx);
	free(pred0); free(oidx); free(ox); free(os); free(odv); free(opred);
	free(ustrat); free(binof); free(bstrat); free(bnum); free(bedlo); free(bedhi); free(xo);
	free(start); free(fill); free(order); free(buf);
	free(xmed); free(predbin);
	free(val); free(obsval); free(obspct); free(simpct);

	info(0, "VPC end\n");
}

/* remove all rows with first column == name, the header line is kept.
 * An empty name deletes the file. */
static void vpc_clean_file(const char* filename, const char* name)
{
	/* no name: delete the whole file */
	if (!name || !*name) {
		remove(filename);
		return;
	}

	/* otherwise rewrite the file without the rows of name into a temp
	 * file, then replace the original with it */
	char tmp[PATH_MAX + 8];
	snprintf(tmp, sizeof(tmp), "%s.tmp", filename);
	var fout = fopen(tmp, "w");
	if (!fout)
		return;

	bool any = false;
	if (!vpc_filter_to(fout, filename, name, &any)) {	/* no file, nothing to clean */
		fclose(fout);
		remove(tmp);
		return;
	}
	if (!vpc_close(fout) || VPC_RENAME(tmp, filename) != 0) {
		remove(tmp);
		fatal(0, "cannot clean VPC file '%s'\n", filename);
	}
}

/// This `pmx_vpc_clean(...)` in openpmxtran function deletes the part
/// of the tables associated with the name. A NULL or empty name deletes
/// the entire VPC output and points files.
void pmx_vpc_clean(const OPENPMX* const pmx, const char* name)
{
	if (!pmx || !pmx->filename || strlen(pmx->filename) == 0)
		return;

	/* clean both the summary file and the points file */
	char filename[PATH_MAX];
	snprintf(filename, sizeof(filename), "%s" OPENPMX_VPC_OUT, pmx->filename);
	vpc_clean_file(filename, name);

	snprintf(filename, sizeof(filename), "%s" OPENPMX_VPC_POINTS, pmx->filename);
	vpc_clean_file(filename, name);
}
