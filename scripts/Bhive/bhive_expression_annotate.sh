#! /bin/bash
### This script makes the figures for the hiv_expression data
### This is specific for the presentation at the end of March 2019

#################### HIV expression vs Absolute Distance to TSS ####################
########## x= Log10(distance to TSS); y= HIV expression; add pearson correlation

##### Turn HIV expression into a bed file
perl /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Bhive/exp_txt2bed.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression.bed

##### Find closest TSS
### Note: get uTSS file from bhive_annotate.sh
module load bedtools/2.26.0
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_uTSS.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS/hiv_expression_closestTSS.bed

### Divide by orientation: same or opposite
awk '$5==$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS/hiv_expression_closestTSS.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS/hiv_expression_closestTSS_same.bed
awk '$5!=$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS/hiv_expression_closestTSS.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS/hiv_expression_closestTSS_opposite.bed

### add newly calculated distance to bed files
awk 'BEGIN {OFS="\t"} {if ($17=="+" && $2>$13) {print $0, $2-$13;} else if($17=="+" && $2>$13) {print $0, $13-$2;} else if($17=="-" && $3>$14) {print $0, $3-$14;} else if($17=="-" && $3<$14) print $0, $14-$3}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS/hiv_expression_closestTSS_same.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS/hiv_expression_closestTSS_same_plus_dis.bed

awk 'BEGIN {OFS="\t"} {if ($17=="+" && $2>$13) {print $0, $2-$13;} else if($17=="+" && $2>$13) {print $0, $13-$2;} else if($17=="-" && $3>$14) {print $0, $3-$14;} else if($17=="-" && $3<$14) print $0, $14-$3}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS/hiv_expression_closestTSS_opposite.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS/hiv_expression_closestTSS_opposite_plus_dis.bed

### Make scatterplot, but manually add the Pearson correlation
module load R/3.3.2-gccmkl
Rscript /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Bhive/BHIVE_expression_scatterplots_closest.R
