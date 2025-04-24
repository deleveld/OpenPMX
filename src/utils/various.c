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

#include <string.h>
#include "c22.h"
#include "various.h"

FILE* results_fopen(const char* name, const char* ext, const char* mode)
{
	let nchars = strlen(name) + strlen(ext) + 1;
	char fname[nchars];
	strcpy(fname, name);
	strcat(fname, ext);
	return fopen(fname, mode);
}

double timespec_time_difference(const struct timespec* const begin,
								const struct timespec* const end)
{
	return (double)(end->tv_sec - begin->tv_sec)*1000. + (double)(end->tv_nsec - begin->tv_nsec)/(double)(1000.*1000.);
}

void timespec_duration(const struct timespec* const begin, double* eval)
{
	struct timespec end;
	clock_gettime(CLOCK_REALTIME, &end);

	const double delta = timespec_time_difference(begin, &end);

	if (*eval <= 0.)
		*eval = delta;
	else
		*eval = (*eval + delta) / 2.;
}

//#define PRINT_MATRIX
#ifdef PRINT_MATRIX
void print_matrix(const gsl_matrix* const p)
{
	if (!p) {
		printf("(null)\n");
		return;
	}
	printf("matrix(%i,%i)\n", (int)p->size1, (int)p->size2);
	forcount(i, p->size1) {
		forcount(j, p->size2)
			printf("%13e ", gsl_matrix_get(p, i, j));
		printf("\n");
	}
}
#endif

//#define PRINT_VECTOR
#ifdef PRINT_VECTOR
void print_vector(const gsl_vector* const p)
{
	if (!p) {
		printf("(null)\n");
		return;
	}
	forcount(i, p->size)
		printf("%13e ", gsl_vector_get(p, i));
	printf("\n");
}
#endif

