# structure of the NONMEM model
NONMEM_MODEL_PREFIX=$(cat <<-MODELPREFIX
	\$PROB Warfarin
	\$INPUT <fromdata>
	\$ABBR DERIV2=NO
MODELPREFIX
)
NONMEM_MODEL_CODE=$(cat <<-MODELCODE
	\$SUBROUTINES ADVAN13 TOL=9
	\$MODEL
		COMP=(ABSORB)
		COMP=(CENTRAL)
		COMP=(TURNOVER)
    \$PK
		size = wt / 70.
		vscale = size**1
		cscale = size**0.75
		kscale = size**(-0.25)
		V = THETA(1) * exp(ETA(1)) * vscale;
		CL = THETA(2) * exp(ETA(2)) * cscale;
		KA = THETA(3) * exp(ETA(3)) * kscale;
		E0 = THETA(5) + ETA(5);
		A_0(3) = E0
		EMAX = THETA(6)
		C50 = THETA(7) * exp(ETA(6))
		KOUT = exp(THETA(8))
	\$DES
		DADT(1) = -A(1)*KA;
		DADT(2) = A(1)*KA - A(2)*(CL/V);
		dcp = A(2)/V;
		deff = dcp/(C50+dcp);
		pd = 1 - EMAX*deff;
		DADT(3) = E0*KOUT*pd - KOUT*A(3)
	\$ERROR
		IPRED = A(2)/V
		BPRED = A(3)
		if (dvid == 1) Y = IPRED*(1+ERR(1)) + ERR(2);
		if (dvid == 2) Y = BPRED + ERR(3);
MODELCODE
)
NONMEM_MODEL_TRUE=$(cat <<-MODELINITIAL
	\$THETA
		(1,   8.04E+00,   20)
		(0,   1.35E-01,   1)
		(0.01,   5.09E-01,   5)
		(0 FIXED)
		(50,   9.65E+01,   150)	; E0
		(0,   1.19E+00,   5)	; EMAX
		(0,   2.19E+00,   5)	; C50
		(-5,   -3.17E+00,   1)	; log(KOUT)
	\$OMEGA 2.61E-02  6.94E-02  4.34E-01 0 FIXED 0 FIXED 0 FIXED
	\$SIGMA 8.91E-03  1.26E+00  2.25E+01
MODELINITIAL
)
NONMEM_MODEL_INITIAL=$(cat <<-MODELINITIAL
	\$THETA
		(1,   8,   20)
		(0,   0.1,   1)
		(0.01,   0.5,   5)
		(0 FIXED)
		(50,   100,   150)	; E0
		(0,   1,   5)	; EMAX
		(0,   2,   5)	; C50
		(-5,   -3,   1)	; log(KOUT)
	\$OMEGA 0.1 0.1 0.1 0 FIXED 0 FIXED 0 FIXED
	\$SIGMA 1 1 100
MODELINITIAL
)

OPENPMX_MODEL_INITIAL=$(cat <<-GRONMEMMODEL
	\$ADVAN(diffeqn_libgsl)
		.nstate = 3,
		.args.diffeqn.abstol = 1e-9,
		.args.diffeqn.reltol = 1e-9,
	\$IMODEL(V, CL, KA, KOUT, E0, EMAX, C50)
		const double size = WT / 70.;
		const double vscale = pow(size, 1);
		const double cscale = pow(size, 0.75);
		const double kscale = pow(size, -0.25);
		V = THETA(1) * exp(ETA(1)) * vscale;
		CL = THETA(2) * exp(ETA(2)) * cscale;
		KA = THETA(3) * exp(ETA(3)) * kscale;
		E0 = THETA(5) + ETA(5);
		A_0(3, E0);
		EMAX = THETA(6);
		C50 = THETA(7) * exp(ETA(6));
		KOUT = exp(THETA(8));
	\$DIFFEQN
		DADT(1) = -A(1)*KA;
		DADT(2) = A(1)*KA - A(2)*(CL/V);
		double dcp = A(2)/V;
		double deff = dcp/(C50+dcp);
		double pd = 1 - EMAX*deff;
		DADT(3) = E0*KOUT*pd - KOUT*A(3);
	\$PREDICT(IPRED, BPRED)
		IPRED = A(2)/V;
		BPRED = A(3);
		if (DVID == 1) Y = IPRED*(1+ERR(1)) + ERR(2);
		if (DVID == 2) Y = BPRED + ERR(3);
	\$THETA
		{ 1,   8,   20, ESTIMATE },
		{ 0,   0.1,   1, ESTIMATE },
		{ 0.01,   0.5,   5, ESTIMATE },
		{ 0, 0, 0, FIXED },
		{ 50,   100,   150, ESTIMATE },	// E0
		{ 0,   1,   5, ESTIMATE },	// EMAX
		{ 0,   2,   5, ESTIMATE },	// C50
		{ -5,   -3,   1, ESTIMATE },	// log(KOUT)
	\$OMEGA(0.1, 0.1, 0.1, 0, 0, 0)
	\$SIGMA(1, 1, 100)
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
	DATASET_HEADER=$(get_datafile_header "models/warfarin_rate0.csv")

	cat >control.${DATASET}.txt <<-CONTROLFILE
	${NONMEM_MODEL_PREFIX}
	\$DATA "models/warfarin_rate0.csv" IGNORE=@
	${NONMEM_MODEL_CODE}
	${NONMEM_MODEL_TRUE}
	\$SIM (${DATASET_SEED}) ONLYSIM
	\$TABLE ${DATASET_HEADER}
	NOPRINT ONEHEADER NOAPPEND file="control.table.txt"
CONTROLFILE
	cat control.${DATASET}.txt
	../utils/do_nonmem_run_singlethread control.${DATASET}.txt
	cat control.${DATASET}.out

	# collect NONMEM results
	gawk '{
		if (NR != 1)
			print
	}' <"control.table.txt" >"simdata/data.${DATASET}.txt"
	rm control.*
	rm gfortran.txt
	rm nmpathlist.txt
}
