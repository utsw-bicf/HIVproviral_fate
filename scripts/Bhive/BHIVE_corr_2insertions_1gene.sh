#! /bin/bash
###This script makes the scatterplots for 2 HIV insertions into 1 gene, comparing expression of insertions

### Get a list of just 2 insertions into 1 gene
Rscript /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Bhive/BHIVE_corr_2insertions_1gene_makelist.R

### filter to put back into R
paste -d " " - - </project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots/delete_sorted_list.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots/delete_sorted_list2.txt

awk -F "\t" '{OFS="\t"}; {if ($17=="+" && $2<$24) {print $15, $11, $33;} else if ($17=="+" && $2>$24) {print $15, $33, $11;} else if ($17=="-" && $2<$24) {print $15, $33, $11;} else if ($17=="-" && $2>$24) print $15, $11, $33;}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots/delete_sorted_list2.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots/two_insertions_1gene.txt

### Put directionality into it (either same, opposite, mixed vs gene)
# same
awk -F "\t" '{OFS="\t"}; {if ($17=="+" && $5=="+" && $27=="+" && $2<$24) {print $15, $11, $33;} else if ($17=="+" && $5=="+" && $27=="+" && $2>$24) {print $15, $33, $11;} else if ($17=="-" && $5=="-" && $27=="-" && $2<$24) {print $15, $33, $11;} else if ($17=="-" && $5=="-" && $27=="-" && $2>$24) print $15, $11, $33;}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots/delete_sorted_list2.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots/two_insertions_1gene_same.txt
# opposite
awk -F "\t" '{OFS="\t"}; {if ($17=="+" && $5=="-" && $27=="-" && $2<$24) {print $15, $11, $33;} else if ($17=="+" && $5=="-" && $27=="-" && $2>$24) {print $15, $33, $11;} else if ($17=="-" && $5=="+" && $27=="+" && $2<$24) {print $15, $33, $11;} else if ($17=="-" && $5=="+" && $27=="+" && $2>$24) print $15, $11, $33;}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots/delete_sorted_list2.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots/two_insertions_1gene_opposite.txt
# mixed
awk -F "\t" '{OFS="\t"}; {if ($17=="+" && $5=="+" && $27=="-" && $2<$24) {print $15, $11, $33;} else if ($17=="+" && $5=="-" && $27=="+" && $2<$24) {print $15, $11, $33;} else if ($17=="+" && $5=="+" && $27=="-" && $2>$24) {print $15, $33, $11;} else if ($17=="+" && $5=="-" && $27=="+" && $2>$24) {print $15, $33, $11;} else if ($17=="-" && $5=="-" && $27=="+" && $2<$24) {print $15, $33, $11;} else if ($17=="-" && $5=="+" && $27=="-" && $2<$24) {print $15, $33, $11;} else if ($17=="-" && $5=="-" && $27=="+" && $2>$24) {print $15, $11, $33;} else if ($17=="-" && $5=="+" && $27=="-" && $2>$24) print $15, $11, $33;}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots/delete_sorted_list2.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots/two_insertions_1gene_mixed.txt

### Then run R to make new scatter plot
Rscript /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Bhive/BHIVE_corr_2insertions_1gene_makescatterplot.R
