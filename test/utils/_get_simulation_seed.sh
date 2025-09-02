set -eu
TEMP=$(mktemp -d)
trap "rm -rf ${TEMP}" EXIT
export LC_ALL=C

SEEDNUM=${1}

SCRIPTNAME=$(basename "${0}")
SCRIPTHOME=$(dirname "${0}")

make_seeds_file()
{
	cp "../../src/utils/c22.h" ${TEMP}

	cat >${TEMP}/seedgen.c <<-ENDSEEDGEN
	#include <stdlib.h>
	#include <stdio.h>
	#include "c22.h"
	#include <gsl/gsl_randist.h>

	int main (void)
	{
		let rng = gsl_rng_alloc(gsl_rng_mt19937);

		/* Birthdays of Joyce, Nette and Isabel :) */
		gsl_rng_set(rng, 200501041406);

		forcount(i, 1000) {	
			let runseed = gsl_rng_uniform_int(rng, 999999);
			printf("%i %lu\n", i + 1, runseed);
		}

		gsl_rng_free(rng);
		return EXIT_SUCCESS;
	}
ENDSEEDGEN
	(cd ${TEMP}; gcc seedgen.c -lgsl -lgslcblas -lm; ./a.out)
}

# make seeds file if not existing yet
if [ ! -f "${SCRIPTNAME}.seeds" ]; then
	make_seeds_file >"${SCRIPTHOME}/${SCRIPTNAME}.seeds"
fi

# get the right seed
gawk '{
	if ($1 == '${SEEDNUM}')
		printf("%lu", $2)
}' <"${SCRIPTHOME}/${SCRIPTNAME}.seeds"

rm -rf ${TEMP}

