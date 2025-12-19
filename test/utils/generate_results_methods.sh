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
	
	MAXNUMBERNODES=$(($(nproc --all) - 4))
	NUMBERNODES="${DO_NONMEM_RUN_NODES:-${MAXNUMBERNODES}}"
	echo nnodes ${NUMBERNODES} >openpmx_nodes.txt

	cat >control.${DATASET}.txt <<-CONTROLFILE
	${NONMEM_MODEL_PREFIX}
	\$DATA "simdata/data.${DATASET}.txt" IGNORE=@
	${NONMEM_MODEL_CODE}
	${NONMEM_MODEL_INITIAL}
	\$ESTM SIG=5 MAX=5000 METHOD=1 INTERACT NOABORT POSTHOC PRINT=1
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
}

###################
# GRONMEM estimation
openpmx()
{
	DATASET=${1}
	RUNNAME=${FUNCNAME[0]}

	MAXNUMBERNODES=$(($(nproc --all) - 4))
	NUMBERNODES="${DO_NONMEM_RUN_NODES:-${MAXNUMBERNODES}}"
	echo nnodes ${NUMBERNODES} >openpmx_nodes.txt

	cat >control.${DATASET}.gr <<-CONTROLFILE
	\$DATA("simdata/data.${DATASET}.txt")
	${OPENPMX_MODEL_INITIAL}
	openpmx.nthread = ${NUMBERNODES};
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

###################
# validate method, same as openpmx
validate()
{
	DATASET=${1}
	RUNNAME=${FUNCNAME[0]}

	MAXNUMBERNODES=$(($(nproc --all) - 4))
	NUMBERNODES="${DO_NONMEM_RUN_NODES:-${MAXNUMBERNODES}}"
	echo nnodes ${NUMBERNODES} >openpmx_nodes.txt

	cat >control.${DATASET}.gr <<-CONTROLFILE
	\$DATA("simdata/data.${DATASET}.txt")
	${OPENPMX_MODEL_INITIAL}
	openpmx.nthread = ${NUMBERNODES};
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

