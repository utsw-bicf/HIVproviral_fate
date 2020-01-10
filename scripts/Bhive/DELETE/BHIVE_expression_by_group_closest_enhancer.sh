#! /bin/bash
#### This script does the nearest enhancer/SuperEnhancer to HIV expression insertion

########## Enhancer
### Group1
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_Enhancer.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same_closestEnhancer.bed

### Group2
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_Enhancer.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite_closestEnhancer.bed

### Group3
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_Enhancer.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite_closestEnhancer.bed

### Group4
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_Enhancer.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene_closestEnhancer.bed

### Group5
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_Enhancer.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene_closestEnhancer.bed

### Group6-8
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6_2gene_opposite_noexp.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group7_2gene_same_noexp.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group8_2gene_opposite_noexp.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6thru8.bed

bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6thru8.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6thru8_2gene.bed

bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6thru8_2gene.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_Enhancer.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6thru8_2gene_closestEnhancer.bed


########## SuperEnhancer
### Group1
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same_closestSuperEnhancer.bed

### Group2
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite_closestSuperEnhancer.bed

### Group3
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite_closestSuperEnhancer.bed

### Group4
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene_closestSuperEnhancer.bed

### Group5
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene_closestSuperEnhancer.bed

### Group6-8
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6thru8_2gene.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6thru8_2gene_closestSuperEnhancer.bed
