#! /bin/sh
### This script makes the metrics_table_sampleprep.xlsx file that contains the NRF, PBC1, PBC2


wdir='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads'
odir='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519'

if [ -e ${odir}/metrics_table_sampleprep.xlsx ]; then
  rm ${odir}/metrics_table_sampleprep.xlsx
fi

files=${wdir}/[GS]*.pbc.qc

echo -e 'Sample\tTotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF\tPBC1\tPBC2' >${odir}/metrics_table_sampleprep.xlsx

for file in ${files}; do
  name=$(basename -z "${file}" .filt.nodup.pbc.qc)
  echo -n ${name} >>${odir}/metrics_table_sampleprep.xlsx
  echo -en '\t' >>${odir}/metrics_table_sampleprep.xlsx
  awk 'FNR>1 {print $0}' ${file} >>${odir}/metrics_table_sampleprep.xlsx
done

