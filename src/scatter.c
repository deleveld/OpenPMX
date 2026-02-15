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
	pthread_cond_t task_is_available;
	pthread_cond_t task_is_done;
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
		pthread_cond_broadcast(&tpool.task_is_available);
		pthread_mutex_unlock(&tpool.mutex);
//		info(0, "tpool: stop ");
		forcount(i, tpool.nrunning) {
			pthread_join(tpool.running[i], 0);
//			info(0, " %i", i);
		}
//		info(0, "\n");
		pthread_mutex_destroy(&tpool.mutex);
		pthread_cond_destroy(&tpool.task_is_available);
		pthread_cond_destroy(&tpool.task_is_done);

		free(tpool.running);

		tpool.init = false;
		tpool.stop = false;
		memset(&tpool, 0, sizeof(tpool));
		memset(&ttasks, 0, sizeof(ttasks));
	}
}

typedef struct {
    const bool is_main_thread;
} THREAD_CONTEXT;

static void* run_threadpool_task(void *_ptr)
{
	var context = (THREAD_CONTEXT*)_ptr;
	bool complete = false;
	while (1) {
		pthread_mutex_lock(&tpool.mutex);

		if (complete) {
			++ttasks.ncomplete;
			/* signal the main thread that a task is done and it can
			 * check if everything is done. */
			pthread_cond_signal(&tpool.task_is_done); 
		}

		/* workers wait here for new individuals */
		while (ttasks.nindivids == 0 && !tpool.stop) {
			if (context && context->is_main_thread)
				goto done;
			pthread_cond_wait(&tpool.task_is_available, &tpool.mutex);
		}
		if (tpool.stop)
			break;

		/* take a task to do, remove from list */
		let individ = ttasks.individs[ttasks.nindivids - 1];
		--ttasks.nindivids;
		pthread_mutex_unlock(&tpool.mutex);

		/* execute task */
		ttasks.threadtask(individ, ttasks.advanfuncs, ttasks.popmodel, 
						  ttasks.nonzero, ttasks.options,
						  ttasks.scatteroptions);
		complete = true;
	}
done:
	pthread_mutex_unlock(&tpool.mutex);
	return _ptr;
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
	let logstream = scatteroptions->logstream;

/// The number of threads is limited to the number of individuals.
/// Otherwise threads will always be waiting and do nothing.
/// The main thread does work as well so we make nthread-1 worker
/// threads.
	var nthread = options->nthread;
	if (nthread > nindivid && scatteroptions->checkout_errors) {
		info(logstream, "nthread (%i) limited to number of individuals (%i)\n", nthread, idata->nindivid);
		nthread = nindivid;
	}
	int nworker = nthread - 1;
	if (nworker <= 0)
		nworker = 0;
#if defined(OPENPMX_PARALLEL_SERIAL)
	nworker = 0;
#endif

#if defined(OPENPMX_PARALLEL_PTHREADS)
	/* if we are already running and the number of threads dont match
	 * then cleanup and init all over again as necessary. */
	if (tpool.init && tpool.nrunning != nworker) 
		scatter_cleanup();
#endif

/// When scattering tasks, we sort individuals based on thier expected
/// runtime. We do the slowest first so we are the most efficient
/// because of lower risk of one long process slowing down finishing the
/// last task. This avoids the "long-tail" inefficiency.
	var index_time = mallocvar(INDIVID_TASK, nindivid);
	forcount(i, nindivid) {
		index_time[i].individ = &individ[i];
		index_time[i].msec = (scatteroptions->stage1_order) ? (individ[i].stage1_msec) : (individ[i].eval_msec);
	}
	qsort(index_time, nindivid, sizeof(INDIVID_TASK), index_time_sort);

	/* individs is allocated here, valid until workers are done */
	var individs = mallocvar(INDIVID*, nindivid);
	forcount(i, nindivid)
		individs[i] = index_time[i].individ;
	free(index_time);

	/* just handle the individuals serially and skip any threading stuff
	 * if we dont use any worker threads. */
	if (nworker == 0) {
		for (int i=0; i<nindivid; i++)
			threadtask(individs[i], advanfuncs, popmodel, nonzero, options, scatteroptions);
		goto done;
	}

#if defined(OPENPMX_PARALLEL_PTHREADS)
	if (!tpool.init) {
		pthread_mutex_init(&tpool.mutex, NULL);
		pthread_cond_init(&tpool.task_is_available, NULL);
		pthread_cond_init(&tpool.task_is_done, NULL);
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
		tpool.running = callocvar(pthread_t, nworker);
		if (!tpool.running)
			fatal(logstream, "tpool: allocation failed");
		info(logstream, "tpool: start (%i)", nthread);

		var ncreated = 0;
		forcount(i, nworker) {
			if (pthread_create(&tpool.running[i], NULL, run_threadpool_task, 0) != 0) {
				warning(logstream, "tpool: create thread failed %d\n", i);
				break;
			}
			info(logstream, " %i", i);
			ncreated++;
		}
		info(logstream, "\n");
		if (ncreated == 0)
			fatal(logstream, "tpool: could not create any threads");
		tpool.nrunning = ncreated;
		tpool.init = true;
	}

	/* start to allow threads to do work, they will lock the mutex to
	 * get task and unlock */
	pthread_cond_broadcast(&tpool.task_is_available);
	pthread_mutex_unlock(&tpool.mutex);
	print_serialize(true);

	/* The main thread must return immediately when there is no work
	 * to be done and wait for the other already started threads to
	 * finish. The other threads will keep running in the background and
	 * wait for new individuals to process. */
	var thread_context = (THREAD_CONTEXT){ .is_main_thread = true };
	run_threadpool_task(&thread_context);

	/* some other threads could still be computing so we wait for them */
	pthread_mutex_lock(&tpool.mutex);
	while (ttasks.ncomplete < nindivid) 
		pthread_cond_wait(&tpool.task_is_done, &tpool.mutex);
	pthread_mutex_unlock(&tpool.mutex);
	assert(ttasks.ncomplete == nindivid);

	print_serialize(false);
	ttasks.individs = 0; /* dont point to invalid memory */
	ttasks.nindivids = 0;

#elif defined(OPENPMX_PARALLEL_OPENMP)
	print_serialize(true);

	/* do this backwards so slowest tasks are done first so we are more
	 * efficient at the end of the individuals */
	omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic, 1)
	for (int i=nindivid-1; i>=0; i--)
		threadtask(individs[i], advanfuncs, popmodel, nonzero, options, scatteroptions);

	print_serialize(false);
#endif

done:
	free(individs);
}

void scatter_cleanup(void)
{
#ifdef OPENPMX_PARALLEL_PTHREADS
	pthreads_cleanup();
#endif
}
