#! /bin/sh

### This bash script makes the files for WashU epigenome browser

### HIV expression
awk 'FNR>1 {OFS="\t"; print $3, $4, $4+1, $1, $10, $5}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.bed

cp /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.bed /archive/shared/DOrso_BICF/WashU_browser_data/
