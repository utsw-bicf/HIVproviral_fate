#! /bin/bash

########## This script takes the overlapping tads and
########## Finds the shortest element for that region
##### This will be done in a "Loop" where you take the tad.bed vs tad.bed,
##### keeping the shortest of A vs B, until you have the smallest file

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad

tad_file='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.tad.2D.bed'

module load bedtools/2.26.0

### Loop 1
bedtools intersect -wao -a ${tad_file} -b ${tad_file} | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round1.bed

awk '{OFS = "\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round1.bed | sort -k 1,1 -k 2,2n | uniq -u >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round1_keep.bed

awk '{OFS = "\t"} {if ($3-$2 < $12-$11) {print $10, $11, $12, $13, $14, $15, $16, $17, $18} else if ($12-$11 < $3-$2) {print $1, $2, $3, $4, $5, $6, $7, $8, $9}}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round1.bed | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round1_filtered.bed


### Loop 2
bedtools intersect -wao -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round1_filtered.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round1_filtered.bed | sort -k 1,1 -k 2,2n | uniq > /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round2.bed

awk '{OFS = "\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round2.bed | sort -k 1,1 -k 2,2n | uniq -u >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round2_keep.bed

awk '{OFS = "\t"} {if ($3-$2 < $12-$11) {print $10, $11, $12, $13, $14, $15, $16, $17, $18} else if ($12-$11 < $3-$2) {print $1, $2, $3, $4, $5, $6, $7, $8, $9}}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round2.bed | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round2_filtered.bed


### Loop 3
bedtools intersect -wao -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round2_filtered.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round2_filtered.bed | sort -k 1,1 -k 2,2n | uniq > /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round3.bed

awk '{OFS = "\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round3.bed | sort -k 1,1 -k 2,2n | uniq -u >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round3_keep.bed

awk '{OFS = "\t"} {if ($3-$2 < $12-$11) {print $10, $11, $12, $13, $14, $15, $16, $17, $18} else if ($12-$11 < $3-$2) {print $1, $2, $3, $4, $5, $6, $7, $8, $9}}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round3.bed | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round3_filtered.bed


######################################################################
######################################################################
######################################################################
########## Ends loops
### Merge results
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/*keep* /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/round3_filtered.bed | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/longest_tad.bed

### Check the results; should be and is the same length as tads
bedtools merge -i /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/longest_tad.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/test_output.bed

########## Average length; 397673
awk '{sum += ($3-$2)} END {print sum/NR}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/longest_tad.bed

### min; 48000
awk 'NR == 1 || ($3-$2) < min {line = $0; min = ($3-$2)}END{print min +1}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/longest_tad.bed

### max; 2013000
awk 'NR == 1 || ($3-$2) > max {line = $0; max = ($3-$2)}END{print max +1}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/longest_tad/longest_tad.bed
