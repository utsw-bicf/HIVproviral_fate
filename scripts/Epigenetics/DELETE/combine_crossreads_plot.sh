#! /bin/sh
###This file merges the *.filt.nodup.taAlign.15.tagAlign.gz.cc.plot.pdf

###Location
wdir="/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/crossReads"

files=${wdir}/[GS]*.plot.pdf

for file in ${files}; do
  outfile=$(basename "${file}" filt.nodup.tagAlign.15.tagAlign.gz.cc.plot.pdf)
  convert -density 300 ${file} ${wdir}/${outfile}.png
done

montage ${wdir}/[GS]*.png -tile 2x13 -mode concatenate ${wdir}/cross_read_plots.pdf

rm ${wdir}/[GS]*.png


