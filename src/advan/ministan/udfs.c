#include "common.h"
#include <math.h>

void calculate_udfs(Config *cfg) {
  /* calculate unit disposition functions */
  int i;
  double prior;
  double temp1, temp2, temp3, temp4;
  double l1, l2, l3, l4;
  int peak_time = 0; /* time when effect site udf peaks              */

  /* Retrieve config variables */
  int delta_seconds = cfg->delta_seconds;
  double *lambda = cfg->lambda;
  double *p_udf = cfg->p_udf;
  double *e_udf = cfg->e_udf;
  double *p_coef = cfg->p_coef;
  double *e_coef = cfg->e_coef;
  /* convert units from /minute to /second */
  double k10 = cfg->k10 / 60.0;
  double k12 = cfg->k12 / 60.0;
  double k13 = cfg->k13 / 60.0;
  double k21 = cfg->k21 / 60.0;
  double k31 = cfg->k31 / 60.0;
  double k41 = cfg->ke0 / 60.0;
  double vc = cfg->vc;

  cube(k10, k12, k21, k13, k31, lambda);

  p_coef[4] = 0;
  lambda[4] = k41;
  if (k31 > 0) {
    p_coef[1] = (k21 - lambda[1]) * (k31 - lambda[1]) /
                (lambda[1] - lambda[2]) / (lambda[1] - lambda[3]) / vc /
                lambda[1];
    p_coef[2] = (k21 - lambda[2]) * (k31 - lambda[2]) /
                (lambda[2] - lambda[1]) / (lambda[2] - lambda[3]) / vc /
                lambda[2];
    p_coef[3] = (k21 - lambda[3]) * (k31 - lambda[3]) /
                (lambda[3] - lambda[2]) / (lambda[3] - lambda[1]) / vc /
                lambda[3];
    e_coef[1] = p_coef[1] / (k41 - lambda[1]) * k41;
    e_coef[2] = p_coef[2] / (k41 - lambda[2]) * k41;
    e_coef[3] = p_coef[3] / (k41 - lambda[3]) * k41;
    e_coef[4] = (k41 - k21) * (k41 - k31) / (lambda[1] - k41) /
                (lambda[2] - k41) / (lambda[3] - k41) / vc;
  } else {
    if (k21 > 0) {
      p_coef[1] = (k21 - lambda[1]) / (lambda[2] - lambda[1]) / vc / lambda[1];
      p_coef[2] = (k21 - lambda[2]) / (lambda[1] - lambda[2]) / vc / lambda[2];
      p_coef[3] = 0;
      e_coef[1] = p_coef[1] / (k41 - lambda[1]) * k41;
      e_coef[2] = p_coef[2] / (k41 - lambda[2]) * k41;
      e_coef[3] = 0;
      e_coef[4] = (k21 - k41) / (lambda[1] - k41) / (lambda[2] - k41) / vc;
    } else {
      p_coef[1] = 1 / lambda[1] / vc;
      p_coef[2] = 0;
      p_coef[3] = 0;
      e_coef[1] = p_coef[1] / (k41 - lambda[1]) * k41;
      e_coef[2] = 0;
      e_coef[3] = 0;
      e_coef[4] = 1 / (lambda[1] - k41) / vc;
    }
  }

  temp1 = 0;
  temp2 = 0;
  temp3 = 0;
  temp4 = 0;

  l1 = exp(-lambda[1]);
  l2 = exp(-lambda[2]);
  l3 = exp(-lambda[3]);
  l4 = exp(-lambda[4]);

  /* calculate udf, plasma concentration, for an infusion of 1/second */
  p_udf[0] = 0;
  for (i = 1; i < 199; i++) {
    temp1 = temp1 * l1 + p_coef[1] * (1 - l1);
    temp2 = temp2 * l2 + p_coef[2] * (1 - l2);
    temp3 = temp3 * l3 + p_coef[3] * (1 - l3);
    p_udf[i] = temp1 + temp2 + temp3;
  }

  /* now calculate udf, effect site, until peak.  Note peak as peak_time */
  temp1 = 0;
  temp2 = 0;
  temp3 = 0;
  temp4 = 0;
  e_udf[0] = 0;

  for (i = 1; i <= delta_seconds; i++) {
    temp1 = temp1 * l1 + e_coef[1] * (1 - l1);
    temp2 = temp2 * l2 + e_coef[2] * (1 - l2);
    temp3 = temp3 * l3 + e_coef[3] * (1 - l3);
    temp4 = temp4 * l4 + e_coef[4] * (1 - l4);
    e_udf[i] = temp1 + temp2 + temp3 + temp4;
  }

  i = delta_seconds;
  prior = e_udf[i - 1];
  while (prior < e_udf[i]) {
    prior = e_udf[i];
    i++;
    temp1 = temp1 * l1;
    temp2 = temp2 * l2;
    temp3 = temp3 * l3;
    temp4 = temp4 * l4;
    e_udf[i] = temp1 + temp2 + temp3 + temp4;
    if (i > 2698) {
      print_error("e_udf i>2698");
      break;
    }
    peak_time = i - 1;
  }
  cfg->peak_time = peak_time;
}
