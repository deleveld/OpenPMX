#ifndef COMMON_H
#define COMMON_H

#include <stddef.h>
#include <stdbool.h>

typedef struct {
	const double k10;
	const double k12;
	const double k21;
	const double k13;
	const double k31;
	const double ke0;
	const double vc;
	const bool target_effect;

	int peak_time;
	double lambda[5];
	double desired;
	double p_coef[5];    /* Virtual coefficients in plasma   */
	double e_coef[5];    /* Virtual coefficients in effect site   */
	double p_udf[200];   /* Udf for plasma, arbitrary infusion length	*/
	double e_udf[2701];  /* Udf for effect site, T = of delta_seconds	*/
	double eff_df[2701]; /* df for effect site, rate = 0 */
	const int delta_seconds;

	/* State accumulators — all zero at t=0 */
	double temp1, temp2, temp3;				/* plasma   */
	double temp1e, temp2e, temp3e, temp4e;	/* effect   */

	/* decay factors */
	double l1, l2, l3, l4;

	double plasma_conc;
	double effect_conc;	
} Config;

void cube(double k10, double k12, double k21, double k13, double k31,
          double *r);
void swap(double *a, double *b);
void calculate_udfs(Config *cfg);
double model(Config *cfg, double temp1, double temp2, double temp3,
             double temp1e, double temp2e, double temp3e, double temp4e,
             double desired);
double virtual_model(Config *cfg, double vm1, double vm2, double vm3,
                     double vm4, int t, int flag);
int find_peak(Config *cfg, int current_time, double rate, double temp1e,
              double temp2e, double temp3e, double temp4e);

void print_double_array(const char *descr, const double *arr, size_t size);
void print_double(const char *descr, const double number);
void print_doubles(const char *descr, const double *arr, size_t size);
void print_error(const char *descr);

#endif
