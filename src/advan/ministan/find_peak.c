#include "common.h"

int find_peak(Config *cfg, int current_time, double rate, double temp1e,
              double temp2e, double temp3e, double temp4e) {
  double current;
  double earlier;
  double later;

  int delta_seconds = cfg->delta_seconds;
  double *e_udf = cfg->e_udf;

  /* set up initial values */
  current =
      virtual_model(cfg, temp1e, temp2e, temp3e, temp4e, current_time, 1) +
      e_udf[current_time] * rate;
  earlier =
      virtual_model(cfg, temp1e, temp2e, temp3e, temp4e, current_time - 1, 1) +
      e_udf[current_time - 1] * rate;
  later =
      virtual_model(cfg, temp1e, temp2e, temp3e, temp4e, current_time + 1, 1) +
      e_udf[current_time + 1] * rate;
  while (current < earlier || current < later) {
    if (current < earlier) {
      if (current_time == delta_seconds) {
        return current_time;
      }
      current_time--;
      later = current;
      current = earlier;
      earlier =
          virtual_model(cfg, temp1e, temp2e, temp3e, temp4e, current_time, 1) +
          e_udf[current_time] * rate;
    } else {
      current_time++;
      earlier = current;
      current = later;
      later = virtual_model(cfg, temp1e, temp2e, temp3e, temp4e,
                            current_time + 1, 1) +
              e_udf[current_time + 1] * rate;
    }
  }
  return current_time;
}