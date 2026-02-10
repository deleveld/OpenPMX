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

/// This file implements a threadpool of waiting processes to process
/// an individual. The task must be of a THREADTASK type and should only
/// touch the data of the individual passed. 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#include "scatter.h"
#include "print.h"
#include "utils/c22.h"

#include "buildflags.h"

/* --------------------------------------------------------------------*/
/// The threadpool is implemented by pthreads or OpenMP. A
/// single-threaded scatter is also possible.
#if defined(OPENPMX_PARALLEL_PTHREADS)
#include <pthread.h>
#elif defined(OPENPMX_PARALLEL_OPENMP)
#include <omp.h>
#endif

typedef struct {
	INDIVID *individ;
	double msec;
} INDIVID_TASK;

static int index_time_sort(const void* a1, const void* a2)
{
	let d1 = (const INDIVID_TASK*)a1;
	let d2 = (const INDIVID_TASK*)a2;
	let t1 = d1->msec;
	let t2 = d2->msec;
	if (t1 < t2)
		return -1;
	if (t1 > t2)
		return 1;
	return 0;
}

#ifdef OPENPMX_PARALLEL_PTHREADS
/* a lot of information taken from
https://nachtimwald.com/2019/04/12/thread-pool-in-c/ */

static struct {
	bool init;
	bool stop;
	pthread_t* running;
	int nrunning;
	pthread_mutex_t mutex;
	pthread_cond_t cond;
} tpool = { };

static struct {
	INDIVID** individs;		/* individuals that need processing */
	int nindivids; 			/* modified by CPUtask as individuals are completed */
	int ncomplete; 			/* count number of completed tasks */
	const ADVANFUNCS* advanfuncs;
	const POPMODEL* popmodel;
	const NONZERO* nonzero;
	const OPTIONS* options;
	const SCATTEROPTIONS* scatteroptions;
	THREADTASK threadtask;
} ttasks;

static void pthreads_cleanup(void)
{
	if (tpool.init) {

		/* signal threads to stop and wait for that */
		pthread_mutex_lock(&tpool.mutex);
		tpool.stop = true;
		pthread_cond_broadcast(&tpool.cond);
		pthread_mutex_unlock(&tpool.mutex);
//		info(logstream, "tpool: stop ");
		forcount(i, tpool.nrunning) {
			pthread_join(tpool.running[i], 0);
//			info(logstream, " %i", i);
		}
//		info(logstream, "\n");
		pthread_mutex_destroy(&tpool.mutex);
		pthread_cond_destroy(&tpool.cond);

		free(tpool.running);

		tpool.init = false;
		tpool.stop = false;
	}
}

static void* run_threadpool_task(void *ptr)
{
	bool complete = false;
	while (1) {
		pthread_mutex_lock(&tpool.mutex);
		if (complete)
			++ttasks.ncomplete;

		/* wait until there is work to do */
		/* main thread (ptr != 0) return when there are no new individuals to get
		 * in that case we still have to wait for the already running threads to finish */
		while (ttasks.nindivids == 0 && !tpool.stop) {
			if (ptr != 0)
				goto done;
			pthread_cond_wait(&tpool.cond, &tpool.mutex);
		}
		if (tpool.stop)
			break;

		/* take the last one in the list, and shorten the list for the other threads */
		let individ = ttasks.individs[ttasks.nindivids - 1];
		--ttasks.nindivids;
		pthread_mutex_unlock(&tpool.mutex);

		/* do the task, we must only touch individual data */
		let popmodel = ttasks.popmodel;
		let nonzero = ttasks.nonzero;
		let options = ttasks.options;
		let advanfuncs = ttasks.advanfuncs;
		let scatteroptions = ttasks.scatteroptions;
		ttasks.threadtask(individ, advanfuncs, popmodel, nonzero, options, scatteroptions);
		complete = true;
	}
done:
    pthread_mutex_unlock(&tpool.mutex);
    return ptr;
}
#endif

/* NOTE: threadtask function must be thread safe and only touch individual data */
void scatter_threads(const IDATA* const idata,
					 const ADVANFUNCS* const advanfuncs,
					 const POPMODEL* const popmodel,
					 const NONZERO* const nonzero,
					 const OPTIONS* const options,
					 SCATTEROPTIONS* scatteroptions,
					 THREADTASK threadtask)
{
	let individ = idata->individ;
	let nindivid = idata->nindivid;

/// When scattering tasks, sort individuals based on thier expected
/// runtime, doing slowest first so we are the most efficient because of
/// lower risk of one long process slowing down finishing the last task.
	var index_time = mallocvar(INDIVID_TASK, nindivid);
	forcount(i, nindivid) {
		index_time[i].individ = &individ[i];
		index_time[i].msec = (scatteroptions->stage1_order) ? (individ[i].stage1_msec) : (individ[i].eval_msec);
	}
	qsort(index_time, nindivid, sizeof(INDIVID_TASK), index_time_sort);
	var individs = mallocvar(INDIVID*, nindivid);
	forcount(i, nindivid)
		individs[i] = index_time[i].individ;
	free(index_time);

/// It does not make sense to start more threads than individuals
/// otherwise threads will always be waiting and do nothing.
	int nthreads = options->nthread;
	if (nthreads < 0)
		nthreads = abs(nthreads);
	if (nthreads > idata->nindivid && scatteroptions->checkout_errors) {
		info(scatteroptions->logstream, "nthread (%i) limited to number of individuals (%i)\n", nthreads, idata->nindivid);
		nthreads = idata->nindivid;
	}

#if defined(OPENPMX_PARALLEL_PTHREADS)
	if (!tpool.init) {
		pthread_mutex_init(&tpool.mutex, NULL);
		pthread_cond_init(&tpool.cond, NULL);
	}

	/* make sure threads wont do anything until we are ready for them */
	pthread_mutex_lock(&tpool.mutex);
	ttasks.individs = individs;
	ttasks.nindivids = nindivid;
	ttasks.ncomplete = 0;
	ttasks.advanfuncs = advanfuncs;
	ttasks.popmodel = popmodel;
	ttasks.nonzero = nonzero;
	ttasks.options = options;
	ttasks.scatteroptions = scatteroptions;
	ttasks.threadtask = threadtask;

	/* startup the thread pool, mutex is locked so we wont start until the
	   condition is signalled */
	if (!tpool.init) {
		tpool.running = callocvar(pthread_t, nthreads);
		info(scatteroptions->logstream, "tpool: start");
		forcount(i, nthreads) {
			pthread_create(&tpool.running[i], NULL, run_threadpool_task, 0);
			info(scatteroptions->logstream, " %i", i);
		}
		info(scatteroptions->logstream, "\n");
		tpool.nrunning = nthreads;
		tpool.init = true;
	}

	/* start to allow threads to do work, they will lock the mutex to get task and unlock */
	pthread_cond_broadcast(&tpool.cond);
	pthread_mutex_unlock(&tpool.mutex);
	print_serialize(true);

/// The main thread does work was well while its waiting.
	/* passing non-zero pointer means we are the main thread, its a bit hacky, I know */
	run_threadpool_task(&ttasks);

	/* some other threads could still be computing so we wait for them */
	bool done = false;
	while (!done) {
		pthread_mutex_lock(&tpool.mutex);
		if (ttasks.ncomplete == nindivid)
			done = true;
		pthread_mutex_unlock(&tpool.mutex);
	}
	assert(ttasks.ncomplete == nindivid);

	print_serialize(false);

#elif defined(OPENPMX_PARALLEL_OPENMP)
	print_serialize(true);

	omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic, 1)
	for (int i=nindivid-1; i>=0; i--)
		threadtask(individs[i], advanfuncs, popmodel, nonzero, options, scatteroptions);

	print_serialize(false);

/* in the absense of any parallel processing we do everything serially */
/* do this forward to be most cache coherent */
#else
	for (int i=0; i<nindivid; i++)
		threadtask(individs[i], advanfuncs, popmodel, nonzero, options, scatteroptions);
#endif

#if defined(OPENPMX_PARALLEL_PTHREADS)
	ttasks.individs = 0; /* dont point to invalid memory */
	ttasks.nindivids = 0;
#endif
	free(individs);
}

void scatter_cleanup(void)
{
#ifdef OPENPMX_PARALLEL_PTHREADS
	pthreads_cleanup();
#endif
}
