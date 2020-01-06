#! /bin/bash
### This script makes the bed files of groups 1-8 to be used for all analysis

### add distance to end of files

########## Group 1: Intergenic, same direction
awk '$5==$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intergenic_closestTSS.bed | sort -k 1,3 -V  | awk '{OFS="\t"}; {if ($5=="+" && $2>$13) {print $0, $2-$13;} else if ($5=="+" && $2<$13) {print $0, $13-$2;} else if ($5=="-" && $3>$14) {print $0, $3-$14;} else if ($5=="-" && $3<$14) print $0, $14-$3}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed

########## Group 2: Intergenic, opposite direction
awk '{if ($17=="+" && $2>$14 && $5!=$17) {print $0;} else if($17=="-" && $2<$13 && $5!=$17) print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intergenic_closestTSS.bed | sort -k 1,3 -V | awk '{OFS="\t"}; {if ($17=="+" && $2>$13) {print $0, $2-$13;} else if ($17=="+" && $2<$13) {print $0, $13-$2;} else if ($17=="-" && $3>$14) {print $0, $3-$14;} else if ($17=="-" && $3<$14) print $0, $14-$3}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed

########## Group 3: Intergenic, opposite direction
awk '{if ($17=="+" && $2<$13 && $5!=$17) {print $0;} else if($17=="-" && $2>$14 && $5!=$17) print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intergenic_closestTSS.bed | sort -k 1,3 -V | awk '{OFS="\t"}; {if ($17=="+" && $2>$13) {print $0, $2-$13;} else if ($17=="+" && $2<$13) {print $0, $13-$2;} else if ($17=="-" && $3>$14) {print $0, $3-$14;} else if ($17=="-" && $3<$14) print $0, $14-$3}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed

########## Group 4: Intragenic, same direction, 1 gene
awk '$5==$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_1gene.bed | sort -k 1,3 -V | awk '{OFS="\t"}; {if ($17=="+" && $2>$13) {print $0, $13-$2, $2-$13;} else if ($17=="+" && $2<$13) {print $0, $13-$2, $13-$2;} else if ($17=="-" && $3>$14) {print $0, $14-$3, $3-$14;} else if ($17=="-" && $3<$14) print $0, $14-$3, $14-$3}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed

########## Group 5: Intragenic, opposite direction, 1 gene
awk '$5!=$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_1gene.bed | sort -k 1,3 -V | awk '{OFS="\t"}; {if ($17=="+" && $2>$13) {print $0, $13-$2, $2-$13;} else if ($17=="+" && $2<$13) {print $0, $13-$2, $13-$2;} else if ($17=="-" && $3>$14) {print $0, $14-$3, $3-$14;} else if ($17=="-" && $3<$14) print $0, $14-$3, $14-$3}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed

########## Group 6: Intragenic, both direction, 2 gene
### Note: this file contains 2 line/insert
awk '{OFS="\t"} {print $1, $2, $3, $5, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes.bed | uniq -u | awk '{print $1, $2, $3, $4}' | sort | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6_2gene_opposite_noexp.bed

########## Group 7: Intragenic, same direction, 2 gene
### Note: this file contains 2 line/insert
comm -23 <(awk '{print $1, $2, $3, $5, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes.bed | sort) <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/delete.bed) | uniq | awk '{OFS="\t"} $4==$5 {print $1, $2, $3, $4}' | sort | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group7_2gene_same_noexp.bed
########## Group 8: Intragenic, opposite direction, 2 gene
### Note: this file contains 2 line/insert
comm -23 <(awk '{print $1, $2, $3, $5, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes.bed | sort) <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/delete.bed) | uniq | awk '{OFS="\t"}$4!=$5 {print $1,$2,$3,$4}' | sort >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group8_2gene_opposite_noexp.bed

### Did not do group 6-8 by hand



#####################################################################################################

