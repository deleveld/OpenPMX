set -eu
TEMP=$(mktemp -d)
trap "rmdir --ignore-fail-on-non-empty $TEMP" EXIT
export LC_ALL=C

gawk '{
	if ($1 == "ITERATION")
		print "DATASET " $0
	else
		print NR-1 $0
}' <"${1}" >"${TEMP}/data"

mv -v "${TEMP}/data" "${1}"
