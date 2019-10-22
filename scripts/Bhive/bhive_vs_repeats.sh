#! /bin/bash

########## This script identifies which HIV insertions are in repeakmasker
HIV='/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression_sorted.bed'

### Download repeatmasker from UCSC genome
awk '{OFS = "\t"} NR >1 {print $6, $7, $8, $12, $11, $10}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/HIVexpression_vs_repeats/human_repeatmasker | sort -k 1,1 -k 2,2n | uniq | perl -p -e 's|\?||g' >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/HIVexpression_vs_repeats/human_repeatmasker.bed

### Merge together
module load bedtools
bedtools intersect -wao -a ${HIV} -b /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/HIVexpression_vs_repeats/human_repeatmasker.bed > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/HIVexpression_vs_repeats/HIV_expression_repeats.bed

### Fix output
awk '{OFS = "\t"} {if ($7 == ".") {print $1, $2, $3, $4, $5, $6, "None", "None", "None", "None", "None", "None", "None"} else {print $0}}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/HIVexpression_vs_repeats/HIV_expression_repeats.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/HIVexpression_vs_repeats/HIV_expression_repeats_fixed.bed

### Draw a pie chart and jitter plot of results in R
