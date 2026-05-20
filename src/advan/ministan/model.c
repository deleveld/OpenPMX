#include "common.h"
#include <math.h>

/* returns rate in units/second */
double model(Config *cfg, double temp1, double temp2, double temp3,
             double temp1e, double temp2e, double temp3e, double temp4e,
             double desired) {
  double rate;    /* calculated pump rate           */
  int temp_peak;  /* temporary peak */
  double current; /* used in effect mode calculation */
  double min_dif; /* minimim acceptable difference */
  double result;
  int i;               /* used in loops */

  /* Retrieve configuration */
  double *eff_df = cfg->eff_df;
  int delta_seconds = cfg->delta_seconds;
  int peak_time = cfg->peak_time;
  double *p_udf = cfg->p_udf;
  double *e_udf = cfg->e_udf;
  int effect_flag = cfg->target_effect;

  int prior_peak_time = peak_time;

  if (effect_flag == 0 ||
      fabs(desired - (temp1e + temp2e + temp3e + temp4e)) < desired * .05)
  /* even in effect mode, if we get very close, we switch and drive the plasma
     again */
  /* this is beneficial both for reasons of processing speed, and to reduce */
  /* the fluctuations seen in pure effect mode */
  {

    /* If first pass, it's easy */
    if (temp1 == 0) {
      return desired / p_udf[delta_seconds];
    } else {
      /* Calculation as described in Bailey/Shafer IEEE */
      result =
          virtual_model(cfg, temp1, temp2, temp3, 0.0, (int)delta_seconds, 0);

      if (desired > result)
        return (desired - result) / p_udf[delta_seconds];
      else
        return 0.0;
    }
  } else {
    /* effect site rate calculation */
    /* first, calculate udf to peak_time */
    /* If first pass, it's easy */
    if (temp1 == 0) {
      prior_peak_time = peak_time;
      return desired / e_udf[peak_time];
    }

    /* zero out old effect site df */
    for (i = 0; i <= peak_time + 1; i++)
      eff_df[i] = 0;

    /* Should the pump be off? */
    if (virtual_model(cfg, temp1e, temp2e, temp3e, temp4e, delta_seconds, 1) >
        desired)
      return 0.0;
    /* minimim acceptable difference */

    min_dif = desired * .0000001; /* simulation should be right on */

    /* Initial settings */
    temp_peak = prior_peak_time;
    if (temp_peak <= delta_seconds)
      temp_peak = delta_seconds + 1;
    rate = (desired -
            virtual_model(cfg, temp1e, temp2e, temp3e, temp4e, temp_peak, 1)) /
           e_udf[temp_peak];
    temp_peak = find_peak(cfg, temp_peak, rate, temp1e, temp2e, temp3e, temp4e);
    current = virtual_model(cfg, temp1e, temp2e, temp3e, temp4e, temp_peak, 1) +
              e_udf[temp_peak] * rate;

    /* Iterate until solution is found */
    while (fabs(current - desired) > min_dif) {
      rate = (desired - virtual_model(cfg, temp1e, temp2e, temp3e, temp4e,
                                      temp_peak, 1)) /
             e_udf[temp_peak];
      temp_peak =
          find_peak(cfg, temp_peak, rate, temp1e, temp2e, temp3e, temp4e);
      current =
          virtual_model(cfg, temp1e, temp2e, temp3e, temp4e, temp_peak, 1) +
          e_udf[temp_peak] * rate;
    }
    prior_peak_time = temp_peak;

    if (rate < .00001 && temp_peak < delta_seconds * 2) {
      rate = 0;
    }
  }
  if (rate > 0)
    return rate;
  else
    return 0.0;
}
