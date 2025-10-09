#
# !!! removed zero observation !!!
# @55,36,0,100,0,Err:502
#

########################################################################
# construct data set
########################################################################

ORIGINAL_DATASET="models/DATASIM.TXT"
DATASET="datasim.data.csv"

# append column to dataset to transform for NONMEM
gawk  'BEGIN{ FS=","; OFS="," };
{
	if (NR == 1) {
		$5 = "RDV"
		printf("%s,DV\n", toupper($0))
	} else {
		if ($5 != 0)
			printf("%s,%f\n", $0, log($5))
	}
}' <"${ORIGINAL_DATASET}" >"${DATASET}"

########################################################################
# model structure methods
########################################################################

# structure of the NONMEM model
NONMEM_MODEL_PREFIX=$(cat <<-MODELPREFIX
	\$PROB Datasim
	\$INPUT <fromdata>
MODELPREFIX
)
NONMEM_MODEL_CODE=$(cat <<-MODELCODE
	\$PRED
		V  = EXP(THETA(1) + ETA(1));
		KE = EXP(THETA(2) + ETA(2));
		KA = KE + EXP(THETA(3) + ETA(3));
		D = 100.;
		CP = (D * KA) / (V * (KA - KE)) * (EXP(-KE * TIME) - EXP(-KA * TIME)) * 100.;
		Y = log(CP) + ERR(1);
MODELCODE
)
NONMEM_MODEL_INITIAL=$(cat <<-MODELINITIAL
	\$THETA
		( 1.61, 3.40, 4.1) ;	/* V */
		( -4.6, -1.2, 0.7) ;	/* KE */
		( -4.6, -0.7, 0.7) ;	/* KA */
	\$OMEGA BLOCK(3)
		1
		0.001 1
		0.001 0.001 1
	\$SIGMA
		1
MODELINITIAL
)

OPENPMX_MODEL_INITIAL=$(cat <<-GRONMEMMODEL
	\$ADVAN(pred)
	\$IMODEL(V, KA, KE)
		V  = exp(THETA(1) + ETA(1));
		KE = exp(THETA(2) + ETA(2));
		KA = KE + exp(THETA(3) + ETA(3));
	\$PREDICT()
		const double dose = 100.;
		const double t = TIME;
		const double CP = (dose * KA) / (V * (KA - KE)) * (exp(-KE * t) - exp(-KA * t)) * 100.;
		Y = log(CP) + ERR(1);
	\$THETA
		{ 1.61, 3.40, 4.1, ESTIMATE },	/* V */
		{ -4.6, -1.2, 0.7, ESTIMATE },	/* KE */
		{ -4.6, -0.7, 0.7, ESTIMATE },	/* KA */
	\$OMEGABLOCK(
		1,
		0.001, 1,
		0.001, 0.001, 1)
	\$SIGMA(1)
	\$MAIN
GRONMEMMODEL
)

###################
# generate dataset
dataset()
{
	DATASET=${1}
	gawk 'BEGIN{ FS="," };
	{
		if (NR == 1)
			print
		if ($1 == '${DATASET}')
			print
	}' <"datasim.data.csv" >"simdata/data.${DATASET}.txt"
}

