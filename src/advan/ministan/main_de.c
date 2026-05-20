#include "common.h"
#include <math.h>
#include <stdio.h>

#include "cube.c"
#include "udfs.c"
#include "helpers.c"
#include "virtual_model.c"
#include "model.c"
#include "find_peak.c"

#include "cfg.c"

int main(int argc, char *argv[]) 
{
	Config cfg = {
		.k10 = 0.0645,	/* Gepts */
		.k12 = 0.1086,
		.k13 = 0.0229,
		.k21 = 0.0245,
		.k31 = 0.0013,
		.ke0 = 0.112,
		.vc  = 14.3,
		.delta_seconds = 10,
		.target_effect = true,
	};

	cfg_init(&cfg);

	printf("step, time_s, rate, plasma_conc, effect_conc\n");

	double desired = 2.0;	/* target concentration (e.g. µg/mL for propofol) */
	int steps = 240;		/* how many delta_seconds intervals to simulate    */
	for (int step = 0; step < steps; step++) {
		if (step == steps/2)
			desired = 3.;

		double rate = calculate_rate(&cfg, desired);
		advance_rate(&cfg, rate);

		int time_s = step * cfg.delta_seconds;
		printf("%d, %d, %.6f, %.6f, %.6f\n",
			   step, time_s, rate, cfg.plasma_conc, cfg.effect_conc);
	}
	return 0;
}

