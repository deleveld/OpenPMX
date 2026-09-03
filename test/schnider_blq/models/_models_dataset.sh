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
		(0.1, 5.85, 25)
		(1, 26, 125)
		(25, 313, 1500)
		(0.1,  1.92, 6)
		(0.1, 1.34, 3)
		(0.1, 0.863, 3)
		(0.01, 0.22, 1)
	\$OMEGA
		1.04E-01
		9.71E-02
		0 FIXED
		1.94E-02
		1.06E-01
		0 FIXED
	\$SIGMA 1 FIXED
	\$ERROR
		IPRED=A(1)/V1
		W = IPRED * THETA(7)
		Y = IPRED + W * ERR(1)
MODELINITIAL
)
NONMEM_MODEL_INITIAL=$(cat <<-MODELINITIAL
	\$THETA
		(0.1, 6, 25)
		(1, 20, 125)
		(25, 200, 800)
		(0.1,  2, 6)
		(0.1, 1, 3)
		(0.1, 1, 3)
		(0.01, 0.22, 1)
	\$OMEGA
		0.1
		0.1
		0 FIXED
		0.1
		0.1
		0 FIXED
	\$SIGMA 1 FIXED
	\$ERROR
		IPRED = A(1)/V1
		W = IPRED * THETA(7)
		Y = IPRED + W * ERR(1)
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
	\$PREDICT(IPRED, W)
		IPRED = A(1)/V1;
		W = IPRED * THETA(7);
		Y = IPRED + W * ERR(1);
	\$THETA
			{ 0.1, 6, 25, ESTIMATE},
			{ 1, 20, 125, ESTIMATE},
			{ 25, 200, 800, ESTIMATE},
			{ 0.1,  2, 6, ESTIMATE},
			{ 0.1, 1, 3, ESTIMATE},
			{ 0.1, 1, 3, ESTIMATE},
			{ 0.01, 0.22, 1, ESTIMATE },
	\$OMEGA(0.1, 0.1, -0, 0.1, 0.1, -0)
	\$SIGMA(-1.)
	\$MAIN
GRONMEMMODEL
)

OPENPMX_BLQ_MODEL_INITIAL=$(cat <<-GRONMEMMODEL
	\$ADVAN(threecomp)
	\$IMODEL(V1, V2, V3, CL, Q2, Q3, SIZE)
		SIZE = WT/70.;
		V1 = THETA(1) * SIZE * exp(ETA(1));
		V2 = THETA(2) * SIZE * exp(ETA(2));
		V3 = THETA(3) * SIZE * exp(ETA(3));
		CL = THETA(4) * pow(SIZE, 0.75) * exp(ETA(4));
		Q2 = THETA(5) * pow(SIZE, 0.75) * exp(ETA(5));
		Q3 = THETA(6) * pow(SIZE, 0.75) * exp(ETA(6));
	\$PREDICT(IPRED, W)
		IPRED = A(1)/V1;
		W = IPRED * THETA(7);
		if (ISBLQ == 0) {
			Y = IPRED + W * ERR(1);
		} else {
			loglik = loglik_llq(BLQVAL, IPRED, W);
		}
	\$THETA
			{ 0.1, 6, 25, ESTIMATE},
			{ 1, 20, 125, ESTIMATE},
			{ 25, 200, 800, ESTIMATE},
			{ 0.1,  2, 6, ESTIMATE},
			{ 0.1, 1, 3, ESTIMATE},
			{ 0.1, 1, 3, ESTIMATE},
			{ 0.01, 0.22, 1, ESTIMATE },
	\$OMEGA(0.1, 0.1, -0, 0.1, 0.1, -0)
	\$SIGMA(-1.)
	\$MAIN
GRONMEMMODEL
)

NONMEM_MODEL_INITIAL_BLQ=$(cat <<-MODELINITIAL
	\$THETA
		(0.1, 6, 25)
		(1, 20, 125)
		(25, 200, 800)
		(0.1,  2, 6)
		(0.1, 1, 3)
		(0.1, 1, 3)
		(0.01, 0.22, 1)
	\$OMEGA
		0.1
		0 FIXED
		0 FIXED
		0 FIXED
		0 FIXED
		0 FIXED
	\$SIGMA 1 FIXED
	\$ERROR
		IPRED = A(1)/V1
		SD = SQRT((IPRED*THETA(7))**2 + 1e-4**2)
		Z = 0
		IF (ISBLQ.EQ.1) THEN
			Z = (BLQVAL-IPRED)/SD
			Y = MAX(PHI(Z), 1e-12)
			F_FLAG = 1
		ELSE
			Y = IPRED + SD*ERR(1)
			F_FLAG = 0
		ENDIF
MODELINITIAL
)

###################
# NONMEM BLQ estimation
nonmem_blq()
{
	DATASET=${1}
	RUNNAME=${FUNCNAME[0]}
	
	MAXNUMBERNODES=$(($(nproc --all) - 4))
	NUMBERNODES="${DO_NONMEM_RUN_NODES:-${MAXNUMBERNODES}}"

	# enable NONMEM M3 method, marking BLQ as observations
R --vanilla --slave <<-RSCRIPT
	d <- read.table("simdata/data.${DATASET}.txt", sep=",", header=TRUE)
	sel <- d\$ISBLQ == 1 & d\$EVID == 2
	d[sel, "EVID"] <- 0
	write.csv(d, file="simdata/data.nonmem_blq.txt", row.names=FALSE, quote=FALSE)
RSCRIPT

	cat >control.${DATASET}.txt <<-CONTROLFILE
	${NONMEM_MODEL_PREFIX}
	\$DATA "simdata/data.nonmem_blq.txt" IGNORE=@
	${NONMEM_MODEL_CODE}
	${NONMEM_MODEL_INITIAL_BLQ}
	\$ESTM SIG=5 MAX=5000 METHOD=1 LAPLACE INTERACT NOABORT POSTHOC PRINT=1
CONTROLFILE
	cat control.${DATASET}.txt
	start=$(date +%s%3N)
	${DO_NONMEM_SCRIPT} "control.${DATASET}.txt" "${NUMBERNODES}"
	end=$(date +%s%3N)
	runtime=$((end - start))
	echo "${DATASET} $runtime" >> "${SCRIPTNAME}.${RUNNAME}.runtime_ms.txt"
    
	# collect NONMEM results and cleanup
	collect_final_estimate "${DATASET}" "control.${DATASET}.ext" "${SCRIPTNAME}.${RUNNAME}.txt"
	
	rm control.*
	rm gfortran.txt
	rm nmpathlist.txt
	rm simdata/data.nonmem_blq.txt
}


###################
# GRONMEM estimation
openpmx_blq()
{
	DATASET=${1}
	RUNNAME=${FUNCNAME[0]}

	MAXNUMBERNODES=$(($(nproc --all) - 4))
	NUMBERNODES="${DO_NONMEM_RUN_NODES:-${MAXNUMBERNODES}}"
	echo nnodes ${NUMBERNODES} >"${RUNNAME}_nodes.txt"

	cat >control.${DATASET}.gr <<-CONTROLFILE
	\$DATA("simdata/data.${DATASET}.txt")
		/* use BLQ as observations */
		if (ISBLQ == 1)
			EVID = 0;
	${OPENPMX_BLQ_MODEL_INITIAL}
//	openpmx.nthread = ${NUMBERNODES};
//	estimate(.nsig=5.);
	estimate();
CONTROLFILE
	cat control.${DATASET}.gr
	start=$(date +%s%3N)
	../../openpmx control.${DATASET}.gr
	end=$(date +%s%3N)
	runtime=$((end - start))
	echo "${DATASET} $runtime" >> "${SCRIPTNAME}.${RUNNAME}.runtime_ms.txt"

	# collect GRONMEM results and cleanup
	collect_final_estimate "${DATASET}" "control.${DATASET}.gr.ext" "${SCRIPTNAME}.${RUNNAME}.txt"
	rm control.*
}

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

	# collect NONMEM results and add BLQ
	Rscript - <<EOF
	df <- read.table("control.table.txt", header=TRUE, skip=1)
	df\$ISBLQ <- 0
	
	# make each occasion a new individual, we dont need EVID=4 anymore
	deltat <- c(TRUE, diff(df\$TIME) < 0)
	newid <- cumsum(deltat)
	df\$ID <- newid
#	sel <- df\$EVID == 4
#	df[sel, "EVID"] <- 1

	# set BLQ value for 20% of dataset
	obsrec <- df[df\$EVID == 0, ]
	blq <- quantile(obsrec[["DV"]], 0.2)

	# update DV
	blqsel <- df\$EVID == 0 & df\$DV <= blq
	df[blqsel, "DV"] <- blq
	df[blqsel, "ISBLQ"] <- 1
	df[blqsel, "EVID"] <- 2
	df\$BLQVAL <- blq

	write.csv(df, "simdata/data.${DATASET}.txt", row.names=FALSE, quote=FALSE)
EOF

	rm control.${DATASET}.out
	rm control.${DATASET}.txt
	rm control.${DATASET}.ext
	rm gfortran.txt
	rm nmpathlist.txt
}
