write_page()
{
	echo "Write: ${1}_${2}"
	gs -sDEVICE=pdfwrite -dFirstPage=${2} -dLastPage=${2} -sOutputFile="${1}_${2}.pdf" -dNumRenderingThreads=8 -dBATCH -dQUIET -dNOPAUSE -dAutoRotatePages=/None ${3}
	gs -sDEVICE=jpeg -r300 -dFirstPage=${2} -dLastPage=${2} -dGraphicsAlphaBits=4 -sOutputFile="${1}_${2}.jpg" -dBATCH -dQUIET -dNOPAUSE -dSAFER ${3}
}

gs -sDEVICE=pdfwrite -sOutputFile="${1}.pdf" -dNumRenderingThreads=8 -dBATCH -dQUIET -dNOPAUSE -dAutoRotatePages=/None ${2}
write_page "${1}.pdf" 1 "${2}"
write_page "${1}.pdf" 2 "${2}"

