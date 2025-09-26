collect_final_estimate()
{
	DATASET=${1}
	EXTFILE=${2}
	OUTFILE=${3}

	# first run needs a header
	NEEDHEADER=0
	if [ ! -f "${OUTFILE}" ]; then
		NEEDHEADER=1
	fi

	# extract the final line
	gawk '{
		if ('${NEEDHEADER}' && $1 == "ITERATION")
			printf("DATASET %s\n", $0)
		if ($1 == -1000000000)
			lastline = $0
	}
	END {
		printf("'${DATASET}' %s\n", lastline)
	}' <"${EXTFILE}" >>"${OUTFILE}"
}

###################
# NONMEM estimation
nonmem()
{
	DATASET=${1}
	RUNNAME=${FUNCNAME[0]}

	cat >control.${DATASET}.txt <<-CONTROLFILE
	${NONMEM_MODEL_PREFIX}
	\$DATA "simdata/data.${DATASET}.txt" IGNORE=@
	${NONMEM_MODEL_CODE}
	${NONMEM_MODEL_INITIAL}
	\$ESTM SIG=6 MAX=5000 METHOD=1 INTERACT NOABORT POSTHOC PRINT=1
CONTROLFILE
	cat control.${DATASET}.txt
	${DO_NONMEM_SCRIPT} control.${DATASET}.txt

	# collect NONMEM results and cleanup
	collect_final_estimate "${DATASET}" "control.${DATASET}.ext" "${SCRIPTNAME}.${RUNNAME}.txt"
	rm control.*
	rm gfortran.txt
	rm nmpathlist.txt
}

###################
# NONMEM estimation (LAPLACE)
nonmemlaplace()
{
	DATASET=${1}
	RUNNAME=${FUNCNAME[0]}

	cat >control.${DATASET}.txt <<-CONTROLFILE
	${NONMEM_MODEL_PREFIX}
	\$DATA "simdata/data.${DATASET}.txt" IGNORE=@
	${NONMEM_MODEL_CODE}
	${NONMEM_MODEL_INITIAL}
	\$ESTM SIG=6 MAX=5000 METHOD=1 INTERACT NOABORT POSTHOC PRINT=1 LAPLACE
CONTROLFILE
	cat control.${DATASET}.txt
	${DO_NONMEM_SCRIPT} control.${DATASET}.txt

	# collect NONMEM results and cleanup
	collect_final_estimate "${DATASET}" "control.${DATASET}.ext" "${SCRIPTNAME}.${RUNNAME}.txt"
	rm control.*
	rm gfortran.txt
	rm nmpathlist.txt
}

###################
# GRONMEM estimation
openpmx()
{
	DATASET=${1}
	RUNNAME=${FUNCNAME[0]}

	cat >control.${DATASET}.gr <<-CONTROLFILE
	\$DATA("simdata/data.${DATASET}.txt")
	${OPENPMX_MODEL_INITIAL}
	estimate();
CONTROLFILE
	cat control.${DATASET}.gr
	../../openpmx control.${DATASET}.gr

	# collect GRONMEM results and cleanup
	collect_final_estimate "${DATASET}" "control.${DATASET}.gr.ext" "${SCRIPTNAME}.${RUNNAME}.txt"
	rm control.*
}

###################
# validate method, same as openpmx
validate()
{
	DATASET=${1}
	RUNNAME=${FUNCNAME[0]}

	cat >control.${DATASET}.gr <<-CONTROLFILE
	\$DATA("simdata/data.${DATASET}.txt")
	${OPENPMX_MODEL_INITIAL}
	estimate();
CONTROLFILE
	cat control.${DATASET}.gr
	../../openpmx control.${DATASET}.gr

	# collect GRONMEM results and cleanup
	collect_final_estimate "${DATASET}" "control.${DATASET}.gr.ext" "${SCRIPTNAME}.${RUNNAME}.txt"
	rm control.*
}

###################
# GRONMEM ICOV estimation
openpmxicov()
{
	DATASET=${1}
	RUNNAME=${FUNCNAME[0]}

	cat >control.${DATASET}.gr <<-CONTROLFILE
	\$DATA("simdata/data.${DATASET}.txt")
	${OPENPMX_MODEL_INITIAL}
	estimate(.stage1.icov_resample = true);
CONTROLFILE
	cat control.${DATASET}.gr
	../../openpmx control.${DATASET}.gr

	# collect GRONMEM results and cleanup
	collect_final_estimate "${DATASET}" "control.${DATASET}.gr.ext" "${SCRIPTNAME}.${RUNNAME}.txt"
	rm control.*
}

###################
# GRONMEM test estimation
openpmxtest()
{
	DATASET=${1}
	RUNNAME=${FUNCNAME[0]}

	cat >control.${DATASET}.gr <<-CONTROLFILE
	\$DATA("simdata/data.${DATASET}.txt")
	${OPENPMX_MODEL_INITIAL}
//	openpmx.brief = true;
	estimate(.stage1.icov_resample=true);
CONTROLFILE
	cat control.${DATASET}.gr
	../../openpmx control.${DATASET}.gr

	# collect GRONMEM results and cleanup
	collect_final_estimate "${DATASET}" "control.${DATASET}.gr.ext" "${SCRIPTNAME}.${RUNNAME}.txt"
	rm control.*
}

###################
# NLMIXR2 
nlmixr2()
{
	DATASET=${1}
	RUNNAME=${FUNCNAME[0]}

	cat >control.${DATASET}.R <<-CONTROLFILE
	df <- read.table("simdata/data.${DATASET}.txt", sep=",", header=TRUE)
	print(head(df))
	${NLMIXR2_MODEL_INITIAL}

	# save EXT like file
	library(readr)
ext_data <- fit\$eta %>%
	mutate(ID = rownames(fit\$eta)) %>%
	relocate(ID)
# Optional: Add individual predicted parameters
	posthoc <- ranef(fit)  # adds individual parameter estimates
	posthoc\$ID <- rownames(posthoc)
	ext_out <- left_join(ext_data, posthoc, by = "ID")
	obj_val <- fit\$objf  # Final objective function value
	ext_out\$OBJ <- obj_val
	ext_out <- ext_out %>% select(ID, OBJ, everything())
	write_csv(ext_out, "control.${DATASET}.R.ext")
	
CONTROLFILE
	cat control.${DATASET}.R
	Rscript control.${DATASET}.R

	dfgdfg

	# collect GRONMEM results and cleanup
	collect_final_estimate "${DATASET}" "control.${DATASET}.R.ext" "${SCRIPTNAME}.${RUNNAME}.txt"
	rm control.*
}

###################
# find if run has been done before
run_exists()
{
	DATASET="${1}"

	if [ -f "${SCRIPTNAME}.${METHODNAME}.txt" ]; then
		gawk '{
			if ($1 == "'${DATASET}'") 
				print "TRUE"
		}' <"${SCRIPTNAME}.${METHODNAME}.txt"
	fi
}

