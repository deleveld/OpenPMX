#include "common.h"
#include <math.h>

void cfg_init(Config* cfg)
{
	calculate_udfs(cfg);

	cfg->l1 = exp(-cfg->lambda[1] * cfg->delta_seconds);
	cfg->l2 = exp(-cfg->lambda[2] * cfg->delta_seconds);
	cfg->l3 = exp(-cfg->lambda[3] * cfg->delta_seconds);
	cfg->l4 = exp(-cfg->lambda[4] * cfg->delta_seconds);
}

double calculate_rate(Config* cfg, const double desired)
{
	for (int i = 0; i <= cfg->peak_time + 1; i++)
		cfg->eff_df[i] = 0;

	return model(cfg,
				 cfg->temp1, cfg->temp2, cfg->temp3,
				 cfg->temp1e, cfg->temp2e, cfg->temp3e, cfg->temp4e,
				 desired);
}

void advance_rate(Config* cfg, const double rate)
{
	cfg->temp1  = cfg->temp1 * cfg->l1 + cfg->p_coef[1] * rate * (1 - cfg->l1);
	cfg->temp2  = cfg->temp2 * cfg->l2 + cfg->p_coef[2] * rate * (1 - cfg->l2);
	cfg->temp3  = cfg->temp3 * cfg->l3 + cfg->p_coef[3] * rate * (1 - cfg->l3);

	cfg->temp1e = cfg->temp1e * cfg->l1 + cfg->e_coef[1] * rate * (1 - cfg->l1);
	cfg->temp2e = cfg->temp2e * cfg->l2 + cfg->e_coef[2] * rate * (1 - cfg->l2);
	cfg->temp3e = cfg->temp3e * cfg->l3 + cfg->e_coef[3] * rate * (1 - cfg->l3);
	cfg->temp4e = cfg->temp4e * cfg->l4 + cfg->e_coef[4] * rate * (1 - cfg->l4);

	cfg->plasma_conc = cfg->temp1 + cfg->temp2 + cfg->temp3;
	cfg->effect_conc = cfg->temp1e + cfg->temp2e + cfg->temp3e + cfg->temp4e;
}

