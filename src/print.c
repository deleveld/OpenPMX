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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>

#include "print.h"
#include "utils/c22.h"

#include "openpmx_compile_options.h"

/* --------------------------------------------------------------------*/
#if defined(OPENPMX_PARALLEL_PTHREADS)
#include <pthread.h>
#elif defined(OPENPMX_PARALLEL_OPENMP)
#include <omp.h>
#endif

/* mutex to syncronize printing from threads */
#if defined(OPENPMX_PARALLEL_PTHREADS)
#define MUTEX_TYPE pthread_mutex_t
#elif defined(OPENPMX_PARALLEL_OPENMP)
#define MUTEX_TYPE omp_lock_t
#elif defined(OPENPMX_PARALLEL_SINGLETHREAD)
#define MUTEX_TYPE int
#endif

static MUTEX_TYPE* thread_print_mutex = 0;

void print_serialize(const bool serial)
{
	/* create and initialize print mutex */
	if (serial) {
		assert(thread_print_mutex == 0);
		MUTEX_TYPE* mutex = mallocvar(MUTEX_TYPE, 1);
		assert(mutex != 0);

#if defined(OPENPMX_PARALLEL_PTHREADS)
		pthread_mutex_init(mutex, NULL);
		pthread_mutex_lock(mutex);
		thread_print_mutex = mutex;
		pthread_mutex_unlock(mutex);
#elif defined(OPENPMX_PARALLEL_OPENMP)
		omp_init_lock(mutex);	/* initial state is unlocked */
		thread_print_mutex = mutex;
#endif

	/* destroy print mutex */
	} else {
		assert(thread_print_mutex != 0);
		MUTEX_TYPE* mutex = thread_print_mutex;

#if defined(OPENPMX_PARALLEL_PTHREADS)
		pthread_mutex_lock(mutex);
		thread_print_mutex = 0;
		pthread_mutex_unlock(mutex);
		pthread_mutex_destroy(mutex);
		
#elif defined(OPENPMX_PARALLEL_OPENMP)
		omp_set_lock(mutex);
		thread_print_mutex = 0;
		omp_unset_lock(mutex);
		omp_destroy_lock(mutex);
#endif
		free(mutex);
	}
}

void openpmx_fputs(FILE* stream1, FILE* stream2, const char* prefix, const char* v)
{
	if (thread_print_mutex) {
#if defined(OPENPMX_PARALLEL_PTHREADS)
		pthread_mutex_lock(thread_print_mutex);
#elif defined(OPENPMX_PARALLEL_OPENMP)
		omp_set_lock(thread_print_mutex);
#endif
	}

	/* actually send the string */
	if (stream1) {
		if (prefix)
			fputs(prefix, stream1);
		fputs(v, stream1);
	}
	if (stream2) {
		if (prefix)
			fputs(prefix, stream2);
		fputs(v, stream2);
	}

	if (thread_print_mutex) {
#if defined(OPENPMX_PARALLEL_PTHREADS)
		pthread_mutex_unlock(thread_print_mutex);
#elif defined(OPENPMX_PARALLEL_OPENMP)
		omp_unset_lock(thread_print_mutex);
#endif
	}
}

static void openpmx_vprintf(FILE* stream1, FILE* stream2, const char* prefix, const char* format, va_list args1)
{
	/* how much space do we need */
    va_list args2;
	va_copy(args2, args1);
    var nchars = vsnprintf(0, 0, format, args2);
    va_end(args2);

	/* reserve some space including zero terminator */
	char* v = mallocvar(char, nchars + 2);
    vsnprintf(v, nchars + 1, format, args1);
	openpmx_fputs(stream1, stream2, prefix, v);
	free(v);
}

void openpmx_printf(FILE* stream1, FILE* stream2, const char* prefix, const char* format, ... )
{
    va_list args1;
    va_start(args1, format);
	openpmx_vprintf(stream1, stream2, prefix, format, args1);
	va_end(args1);
}

static const char openpmx_error_prefix[] = "error: ";
static const char openpmx_warning_prefix[] = "warning: ";
static const char openpmx_info_prefix[] = "";

void fatal(FILE* stream, const char* format, ... )
{
    va_list args1;
    va_start(args1, format);
	openpmx_vprintf(stream, stderr, openpmx_error_prefix, format, args1);
	va_end(args1);
	exit(EXIT_FAILURE);
}

void warning(FILE* stream, const char* format, ... )
{
    va_list args1;
    va_start(args1, format);
	openpmx_vprintf(stream, stderr, openpmx_warning_prefix, format, args1);
	va_end(args1);
}

void info(FILE* stream, const char* format, ... )
{
    va_list args1;
    va_start(args1, format);
	openpmx_vprintf(stream, stdout, openpmx_info_prefix, format, args1);
	va_end(args1);
}

