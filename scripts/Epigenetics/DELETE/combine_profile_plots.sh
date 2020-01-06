#! /bin/sh
### This script combines multiple profile plots onto 1 image

###Location
wdir="/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/plotProfile"

files=${wdir}/[GS]*.pdf

for file in ${files}; do
  outfile=$(basename "${file}" .pdf)
  convert -density 300 ${file} ${wdir}/${outfile}.png
done

montage ${wdir}/[GS]*.png -tile 2x13 -mode concatenate ${wdir}/profile_plots.pdf

rm ${wdir}/[GS]*.png
