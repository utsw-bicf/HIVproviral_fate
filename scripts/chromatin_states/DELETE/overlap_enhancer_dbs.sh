#! /bin/bash

########## This script finds the enhancers and superenhancers, from rose and homer
########## and finds overlaps (end-to-end or in-both)

module load bedtools

### enhancers
# Homer: 11086
# Rose: 292
# O/L: 208
he='/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM/H_enhancers.bed'
re='/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM/R_enhancers.bed'
Odir='/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM'

# end-to-end
bedtools intersect -wo -a ${he} -b ${re} | awk '{OFS="\t"}; {if ($2 < $5) {start = $2} else {start = $5}; if ($3 > $6) {end=$3} else{end=$6};  {print $1,start,end}}' >${Odir}/enhancers_e2e.txt

# only overlapping regions
bedtools intersect -a ${he} -b ${re} >${Odir}/enhancers_ol.txt

### super enhancers
# Homer: 736
# Rose: 4792
# ENCODE: 929
# O/L H&R: 369
# O/L all3 e2e: 106
# O/L all3 ol: 21
hse='/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM/H_superenhancers.bed'
rse='/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM/R_superenhancers.bed'
ese='/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM/ENCODE_superenhancers.bed'
Odir='/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM'

# end-to-end
bedtools intersect -wo -a ${hse} -b ${rse} | awk '{OFS="\t"}; {if ($2 < $5) {start = $2} else {start = $5}; if ($3 > $6) {end=$3} else{end=$6};  {print $1,start,end}}' >${Odir}/superenhancers_e2e.txt
bedtools intersect -wo -a ${Odir}/superenhancers_e2e.txt -b ${ese} | awk '{OFS="\t"}; {if ($2 < $5) {start = $2} else {start = $5}; if ($3 > $6) {end=$3} else{end=$6};  {print $1,start,end}}' >${Odir}/superenhancers_e2e_all3.txt


# only overlapping regions
bedtools intersect -a ${he} -b ${re} >${Odir}/superenhancers_ol.txt

bedtools intersect -a ${Odir}/superenhancers_ol.txt -b ${ese} >${Odir}/superenhancers_ol_all3.txt

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
### add to chromHMM histone only, reorder overlapping
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder/combine_db
