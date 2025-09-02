# structure of the NONMEM model
NONMEM_MODEL_PREFIX=$(cat <<-MODELPREFIX
	\$PROB Fentanyl
	\$INPUT <fromdata>
MODELPREFIX
)
NONMEM_MODEL_CODE=$(cat <<-MODELCODE
	\$SUBROUTINES ADVAN11 TRANS4
	\$PK
		TH1 = THETA(1)
		TH2 = THETA(2)
		TH3 = THETA(3)
		TH4 = THETA(4)
		TH5 = THETA(5)
		TH6 = THETA(6)
		TH7 = THETA(7)
		TH8 = THETA(8)
		V1 = TH1*(WT/70)**TH7*EXP(ETA(1)) 
		V2 = TH2*(WT/70)**TH7*EXP(ETA(2))
		V3 = TH3*(WT/70)**TH7*EXP(ETA(3)) 
		CL = TH4*(WT/70)**TH8*EXP(ETA(4))
		Q2 = TH5*(WT/70)**TH8*EXP(ETA(5))
		Q3 = TH6*(WT/70)**TH8*EXP(ETA(6)) 
		S1=V1
MODELCODE
)
# TODO: Change this to the exact estimated values, these come from
# a refit, but we should use the ones from the publication
NONMEM_MODEL_TRUE=$(cat <<-MODELINITIAL
	\$THETA
		(0, 1.01E+01, 25)
		(1, 2.63E+01, 125)
		(25, 2.05E+02, 500)
		(0, 7.03E-01, 2)
		(0, 2.36E+00, 5)
		(0, 1.49E+00, 3)
		(0, 1.23E+00, 3)
		(0, 3.10E-01, 1)
	\$OMEGA BLOCK(6)
		4.57E-01
		2.80E-01  5.16E-01
		8.84E-02  1.15E-01  1.35E-01
		4.68E-02  1.31E-02  7.97E-03  9.46E-02
		2.00E-01  3.43E-01  5.89E-02  2.89E-02  4.51E-01
		8.26E-02 -5.14E-02 -1.91E-02  2.26E-02  1.61E-01  2.45E-01
	\$SIGMA 2.62E-02
	\$ERROR
		IPRED=A(1)/V1
		Y = IPRED*(1 + ERR(1))
MODELINITIAL
)
NONMEM_MODEL_INITIAL=$(cat <<-MODELINITIAL
	\$THETA
		(0, 10, 25)
		(1, 25, 125)
		(25, 200, 500)
		(0, 0.7, 2)
		(0, 2, 5)
		(0, 1.4, 3)
		(0, 1, 3)
		(0, 0.75, 1)
	\$OMEGA BLOCK(6)
		.5
		0.001 .5
		0.001 0.001 .5
		0.001 0.001 0.001 .2
		0.001 0.001 0.001 0.001 .2
		0.001 0.001 0.001 0.001 0.001 .2
	\$SIGMA 0.2
	\$ERROR
		IPRED=A(1)/V1
		Y = IPRED*(1 + ERR(1));
MODELINITIAL
)

OPENPMX_MODEL_INITIAL=$(cat <<-GRONMEMMODEL
	\$ADVAN(threecomp)
	\$IMODEL(V1, V2, V3, CL, Q2, Q3, SIZE)
		const double TH1 = THETA(1);
		const double TH2 = THETA(2);
		const double TH3 = THETA(3);
		const double TH4 = THETA(4);
		const double TH5 = THETA(5);
		const double TH6 = THETA(6);
		const double TH7 = THETA(7);
		const double TH8 = THETA(8);
		V1 = TH1*pow(WT/70.,TH7)*exp(ETA(1));
		V2 = TH2*pow(WT/70.,TH7)*exp(ETA(2));
		V3 = TH3*pow(WT/70.,TH7)*exp(ETA(3));
		CL = TH4*pow(WT/70.,TH8)*exp(ETA(4));
		Q2 = TH5*pow(WT/70.,TH8)*exp(ETA(5));
		Q3 = TH6*pow(WT/70.,TH8)*exp(ETA(6));
	\$PREDICT(IPRED)
		IPRED = A(1)/V1;
		Y = IPRED*(1 + ERR(1));
	\$THETA
		{ 0, 10, 25, ESTIMATE },
		{ 1, 25, 125, ESTIMATE },
		{ 25, 200, 500, ESTIMATE },
		{ 0, 0.7, 2, ESTIMATE },
		{ 0, 2, 5, ESTIMATE },
		{ 0, 1.4, 3, ESTIMATE },
		{ 0, 1, 3, ESTIMATE },
		{ 0, 0.75, 1, ESTIMATE },
	\$OMEGABLOCK(
		.5,
		0.001, .5,
		0.001, 0.001, .5,
		0.001, 0.001, 0.001, .2,
		0.001, 0.001, 0.001, 0.001, .2,
		0.001, 0.001, 0.001, 0.001, 0.001, .2)
	\$SIGMA(.2)
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
	DATASET_HEADER=$(get_datafile_header "./models/bae_fentanyl.csv")

	cat >control.${DATASET}.txt <<-CONTROLFILE
	${NONMEM_MODEL_PREFIX}
	\$DATA "models/bae_fentanyl.csv" IGNORE=@
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
