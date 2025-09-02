construct_blank_datafile()
{
R --vanilla --slave <<-RSCRIPT
	if (!file.exists("simdata/${SCRIPTNAME}.data.csv")) {
		d <- data.frame(ID=integer(),
						TIME=numeric(),
						DV=numeric(),
						AMT=numeric(),
						EVID=numeric())

		make_individual_data <- function(id, dose)
		{
			recordtime <- c(0, sort(round(sample(runif(600, 0, 1440), 4) / 60, 2)))

			individ <- d[0, ]
			individ[1:length(recordtime),"ID"] <- id
			individ\$TIME <- recordtime
			individ\$DV <- ifelse(recordtime == 0, 0, 1)
			individ\$AMT <- ifelse(recordtime == 0, dose, 0)
			individ\$EVID <- ifelse(recordtime == 0, 1, 0)

			individ
		}

		filename <- sprintf("simdata/${SCRIPTNAME}.data.csv")
		cat(sprintf("... write %s\n", filename))
		write.table(d, file=filename, row.names=FALSE, col.names=TRUE, quote=FALSE)

		id <- 1
		for (idnum in 1:30) {
			idata <- make_individual_data(id, 10)
			write.table(idata, file=filename, row.names=FALSE, col.names=FALSE, append=TRUE)
			id <- id + 1
		}
		for (idnum in 1:30) {
			idata <- make_individual_data(id, 30)
			write.table(idata, file=filename, row.names=FALSE, col.names=FALSE, append=TRUE)
			id <- id + 1
		}
		for (idnum in 1:30) {
			idata <- make_individual_data(id, 60)
			write.table(idata, file=filename, row.names=FALSE, col.names=FALSE, append=TRUE)
			id <- id + 1
		}
		for (idnum in 1:30) {
			idata <- make_individual_data(id, 120)
			write.table(idata, file=filename, row.names=FALSE, col.names=FALSE, append=TRUE)
			id <- id + 1
		}
	}
RSCRIPT
}

# structure of the NONMEM model
NONMEM_MODEL_PREFIX=$(cat <<-MODELPREFIX
	\$PROB Sparse
	\$INPUT <fromdata>
	\$ABBR DERIV2=NO
MODELPREFIX
)
NONMEM_MODEL_CODE=$(cat <<-MODELCODE
	\$SUBROUTINES ADVAN2 TRANS2
	\$PK
		CL = EXP(THETA(1)+ETA(1))
		V  = EXP(THETA(2)+ETA(2))
		KA = EXP(THETA(3)+ETA(3))
		S2 = V
	\$ERROR
		IPRED = A(2)/V
		Y     = IPRED*(1+ERR(1))
MODELCODE
)
NONMEM_MODEL_TRUE=$(cat <<-MODELINITIAL
	\$THETA	(1.38629436112 FIXED)		;CL
			(4.24849524205 FIXED)		;V
			(0 FIXED)	;Ka
	\$OMEGA   0.09 0.09 0.09
	\$SIGMA   0.04
MODELINITIAL
)
NONMEM_MODEL_INITIAL=$(cat <<-MODELINITIAL
	\$THETA	(-1, 2, 4)		;CL
			(2, 4, 6)		;V
			(-0.5, -0.1, 0.5)	;Ka
	\$OMEGA 0.2 0.2 0.2
	\$SIGMA 0.1
MODELINITIAL
)

OPENPMX_MODEL_INITIAL=$(cat <<-GRONMEMMODEL
	\$ADVAN(onecomp_depot)
	\$IMODEL(V, CL, KA, RESSD)
		CL = exp(THETA(1) + ETA(1));
		V  = exp(THETA(2) + ETA(2));
		KA = exp(THETA(3) + ETA(3));
	\$PREDICT(IPRED)
		IPRED = A(2)/V;
		Y = IPRED * (1 + ERR(1));
	\$THETA
		{   -1,   2,   4,  ESTIMATE     },	/* true CL=log(4) */
		{   2,    4,   6,  ESTIMATE     },	/* true V =log(70) */
		{   -0.5,  -0.1,   0.5,  ESTIMATE     },	/* true KA=log(1) */
	\$OMEGA(0.2, 0.2, 0.2)
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
	construct_blank_datafile

	DATASET=${1}

	DATASET_SEED=$(../utils/_get_simulation_seed.sh ${DATASET})
	DATASET_HEADER=$(get_datafile_header "simdata/${SCRIPTNAME}.data.csv")

	cat >control.${DATASET}.txt <<-CONTROLFILE
	${NONMEM_MODEL_PREFIX}
	\$DATA "simdata/${SCRIPTNAME}.data.csv" IGNORE=@
	${NONMEM_MODEL_CODE}
	${NONMEM_MODEL_TRUE}
	\$SIM (${DATASET_SEED}) ONLYSIM
	\$TABLE ${DATASET_HEADER} NOPRINT ONEHEADER NOAPPEND file="control.table.txt"
CONTROLFILE
	cat control.${DATASET}.txt
	../utils/do_nonmem_run control.${DATASET}.txt

	# collect NONMEM results
	gawk '{
		if (NR != 1)
			print
	}' <"control.table.txt" >"simdata/data.${DATASET}.txt"
	rm control.*
	rm gfortran.txt
	rm nmpathlist.txt
}

