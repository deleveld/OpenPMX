NONMEM_MODEL_PREFIX=$(cat <<-MODELPREFIX
	\$PROB Schnider
	\$INPUT <fromdata>
MODELPREFIX
)
NONMEM_MODEL_CODE=$(cat <<-MODELCODE
	\$SUBROUTINES ADVAN11 TRANS4
	\$PK
		SIZE = WT/70.;
		V1 = THETA(1) * SIZE * exp(ETA(1));
		V2 = THETA(2) * SIZE * exp(ETA(2));
		V3 = THETA(3) * SIZE * exp(ETA(3));
		CL = THETA(4) * (SIZE**0.75) * exp(ETA(4));
		Q2 = THETA(5) * (SIZE**0.75) * exp(ETA(5));
		Q3 = THETA(6) * (SIZE**0.75) * exp(ETA(6));
		S1=V1
MODELCODE
)
NONMEM_MODEL_TRUE=$(cat <<-MODELINITIAL
	\$THETA
		(0, 5.85, 25)
		(1, 26, 125)
		(25, 313, 500)
		(0,  1.92, 6)
		(0, 1.34, 3)
		(0, 0.863, 3)
	\$OMEGA
		1.04E-01
		9.71E-02
		0 FIXED
		1.94E-02
		1.06E-01
		0 FIXED
	\$SIGMA 5.15E-02
	\$ERROR
		IPRED=A(1)/V1
		Y = IPRED*(1 + ERR(1));
MODELINITIAL
)
NONMEM_MODEL_INITIAL=$(cat <<-MODELINITIAL
	\$THETA
		(0, 6, 25)
		(1, 20, 60)
		(50, 200, 500)
		(0.1, 2, 6)
		(0.1, 1, 4)
		(0.1, 0.5, 3)
	\$OMEGA
		0.1
		0.1
		0 FIXED
		0.1
		0.1
		0 FIXED
	\$SIGMA 0.1
	\$ERROR
		IPRED=A(1)/V1
		Y = IPRED*(1 + ERR(1));
MODELINITIAL
)

OPENPMX_MODEL_INITIAL=$(cat <<-GRONMEMMODEL
	\$ADVAN(threecomp)
	\$IMODEL(V1, V2, V3, CL, Q2, Q3, SIZE)
		SIZE = WT/70.;
		V1 = THETA(1) * SIZE * exp(ETA(1));
		V2 = THETA(2) * SIZE * exp(ETA(2));
		V3 = THETA(3) * SIZE * exp(ETA(3));
		CL = THETA(4) * pow(SIZE, 0.75) * exp(ETA(4));
		Q2 = THETA(5) * pow(SIZE, 0.75) * exp(ETA(5));
		Q3 = THETA(6) * pow(SIZE, 0.75) * exp(ETA(6));
	\$PREDICT(IPRED)
		IPRED = A(1)/V1;
		Y = IPRED*(1 + ERR(1));
	\$THETA
			{ 0, 6, 25, ESTIMATE},
			{ 1, 20, 125, ESTIMATE},
			{ 25, 200, 500, ESTIMATE},
			{ 0,  2, 6, ESTIMATE},
			{ 0, 1, 3, ESTIMATE},
			{ 0, 1, 3, ESTIMATE},
	\$OMEGA(0.1, 0.1, -0, 0.1, 0.1, -0)
	\$SIGMA(0.1)
	\$MAIN
GRONMEMMODEL
)

get_datafile_header() {
	gawk '{
		printf("%s", toupper($0))
		exit 0
	}' ${1}
}

###################
# generate dataset
dataset()
{
	DATASET=${1}

	DATASET_SEED=$(../utils/_get_simulation_seed.sh ${DATASET})
	DATASET_HEADER=$(get_datafile_header "models/data.csv")

	cat >control.${DATASET}.txt <<-CONTROLFILE
	${NONMEM_MODEL_PREFIX}
	\$DATA "models/data.csv" IGNORE=@
	${NONMEM_MODEL_CODE}
	${NONMEM_MODEL_TRUE}
	\$SIM (${DATASET_SEED}) ONLYSIM
	\$TABLE ${DATASET_HEADER}
	NOPRINT ONEHEADER NOAPPEND file="control.table.txt"
CONTROLFILE
	cat control.${DATASET}.txt
	../utils/do_nonmem_run control.${DATASET}.txt
	cat control.${DATASET}.out

	# collect NONMEM results
	gawk '{
		if (NR != 1)
			print
	}' <"control.table.txt" >"simdata/data.${DATASET}.txt"
	rm control.${DATASET}.out
	rm control.${DATASET}.txt
	rm control.${DATASET}.ext
	rm gfortran.txt
	rm nmpathlist.txt
}
