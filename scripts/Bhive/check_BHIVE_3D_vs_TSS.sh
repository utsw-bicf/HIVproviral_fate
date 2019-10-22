#! /bin/bash

########## This script checks to make sure that nearest TSS is
########## still inside the same loop as the HIV insertion/expression
########## This also checks to see if still inside the same laminB sub-compartment

of='/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots/check_3D_of_TSS'


tad='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/shortest_tad/shortest_tad.bed'
lamin='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments_hg38.bed'


g1='/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed'
g2='/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed'
g3='/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed'
g4='/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed'
g5='/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed'

cat ${g1} ${g2} ${g3} ${g4} ${g5} >${of}/HIVexpression_nearestTSS.bed

### Check against tads

## Fix tad coordinates
awk '{OFS = "\t"} {if ($3<$2) {print $1, $2, $3, $4, $5, $6} else {print $0}}' ${tad} >${of}/shortest_tad.bed


module load bedtools/2.26.0
bedtools intersect -wao -a ${of}/HIVexpression_nearestTSS.bed -b ${of}/shortest_tad.bed >${of}/HIVexpression_nearestTSS_tad.bed

awk -F "\t" '{OFS = "\t"} {if ($24 ~ /chr/ && $17 == "+" && ($13 < $25 || $13 > $26)) {print $0} else if ($24 ~ /chr/ && $17 == "+" && ($14 < $25 || $14 > $26)) {print $0}}' ${of}/HIVexpression_nearestTSS_tad.bed >${of}/HIVexpression_nearestTSS_tad_badresults.bed


### Check against lamin
bedtools intersect -wao -a ${of}/HIVexpression_nearestTSS.bed -b ${lamin} >${of}/HIVexpression_nearestTSS_lamin.bed

awk -F "\t" '{OFS = "\t"} {if ($24 ~ /chr/ && $17 == "+" && ($13 < $25 || $13 > $26)) {print $0} else if ($24 ~ /chr/ && $17 == "+" && ($14 < $25 || $14 > $26)) {print $0}}' ${of}/HIVexpression_nearestTSS_lamin.bed >${of}/HIVexpression_nearestTSS_lamin_badresults.bed


### How many overlap between badresults Tad and lamin
awk '{OFS = "\t"} {print $1, $2, $3}' ${of}/HIVexpression_nearestTSS_tad_badresults.bed | sort -k 1,1 -k 2,2n >${of}/HIVexpression_nearestTSS_tad_badresults_3.bed
awk '{OFS = "\t"} {print $1, $2, $3}' ${of}/HIVexpression_nearestTSS_lamin_badresults.bed | sort -k 1,1 -k 2,2n >${of}/HIVexpression_nearestTSS_lamin_badresults_3.bed

bedtools intersect -wa -a ${of}/HIVexpression_nearestTSS_tad_badresults_3.bed -b ${of}/HIVexpression_nearestTSS_lamin_badresults_3.bed >${of}/badresults_inboth.bed
