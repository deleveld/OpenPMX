#include "common.h"
#include <math.h>

double virtual_model(Config *cfg, double vm1, double vm2, double vm3,
                     double vm4, int t, int flag) {
  double temp;
  double vmf1, vmf2, vmf3, vmf4;

  int peak_time = cfg->peak_time;
  double *eff_df = cfg->eff_df;
  double *lambda = cfg->lambda;

  if (flag && t > peak_time + 1) {
    print_error("Time error in virtual_model");
    return 0.0;
  }

  if (flag) {
    if (eff_df[t] > 0) {
      return eff_df[t];
    }
  }
  if ((lambda[1] * t) > 100.0)
    vmf1 = 0;
  else
    vmf1 = exp(-lambda[1] * t);

  if ((lambda[2] * t) > 100.0)
    vmf2 = 0;
  else
    vmf2 = exp(-lambda[2] * t);

  if ((lambda[3] * t) > 100.0)
    vmf3 = 0;
  else
    vmf3 = exp(-lambda[3] * t);

  if ((lambda[4] * t) > 100.0)
    vmf4 = 0;
  else
    vmf4 = exp(-lambda[4] * t);

  temp = vm1 * vmf1 + vm2 * vmf2 + vm3 * vmf3 + vm4 * vmf4;

  if (flag)
    eff_df[t] = temp;
  return temp;
}