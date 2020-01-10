#! /bin/bash

########## This script calculates the geometric mean of features
########## Features include: enhancers, protein-coding genes, 7 histones

module load bedtools
module load samtools/gcc/1.6

##### Enhancers; 2180 total
### Get bed file of enhancers
cat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/H3K27ac_enhancers.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/H3K4me3_enhancers.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/H3K4me1_enhancers.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/TTseq_enhancers.bed | sort -k 1,1 -k 2,2n | bedtools merge -i stdin >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/enhancers.bed

### Merge bam files
samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_merged.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.dedup.bam
samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_merged.bam

### Try feature counts instead
# make bed file SAF file
awk '{OFS = "\t"} {print "enhancer" NR, $1, $2, $3, "."}' /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/enhancers.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/enhancers.saf

featureCounts --ignoreDup \
  -T 16 \
  -p \
  -F SAF \
  -a /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/enhancers.saf \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/TTseq_merged_FC.cts \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_merged.bam

# Reorder file
awk '{OFS = "\t"} NR >1 {print $2, $3, $4, $1, $6, $7}' /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/TTseq_merged_FC.cts | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/TTseq_merged_FC.bed

### Find HIV in tad
bedtools intersect -wo -a /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_TadLoop/Tad.bed | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/HIVexpression_inTad.bed

### Calculate in tad
#cat /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/enhancer_head1.bed | while read line; do
cat input.bed | while read line; do
  echo $line | perl -p -e 's| |\t|g' | awk '{OFS = "\t"} {print $7, $8, $9, $1, $2, $3, $4, $5, $6}' >delete1.bed
  #awk '{OFS = "\t"} {print $7, $8, $9, $1, $2, $3, $4, $5, $6}' ${line} >delete1.bed
  bedtools intersect -wao -a delete1.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/TTseq_merged_FC.bed >delete2.bed
done

chr1	9593999	10493998	chr1	10096116	10096117	ATTGTTATACTTTTTCCGTC	0.69142326908866	+
chr11	83133800	83134000	chr11	82679999	83828998	200
chr11	83133800	83134000	chr11	82862999	83291998	200




######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
### Don't use

### Get coverage of each base in merged bam file; only for enhancer regions
# Format enhancers

#module load gatk/3.8
#java -jar /cm/shared/apps/gatk/3.7/GenomeAnalysisTK.jar \
#   -T DepthOfCoverage \
#   -R /project/shared/bicf_workflow_ref/human/GRCh38/genome.fa \
#   -o /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/TTseq_merged_GATKcounts.bed \
#   -I /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_merged.bam \
#   -L /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/enhancers.bed \
#   -U ALLOW_N_CIGAR_READS

### Get coverage of each base in merged bam file
#bedtools genomecov -d -ibam /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_merged.bam -g /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/TTseq_merged_counts.bed
# Make 4 column
#awk '{OFS = "\t"} {print $1, $2-1, $3, $3}' /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/TTseq_merged_counts.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/TTseq_merged_counts_fixed.bed

### Remove zero values; remove alternate chromosomes
#awk '$3 > 0 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/TTseq_merged_counts.bed | awk '$1 ~ /chr/ {print $0}' >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/TTseq_merged_nonzero_counts.bed

### Go line-by-line in bed file and calculate geometric mean
#cat /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/enhancer_head1.bed | while read line; do
#  echo $line | perl -p -e 's| |\t|g' >input.bed
  # Get overlaps
#  bedtools intersect -wao -a input.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/TTseq_merged_counts.bed >test.bed
#done

### Go line-by-line in bed file and calculate geometric mean
#perl /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/machine_learning/calculate_geometric_mean.pl /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/TTseq_merged_counts.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/enhancers.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/Geometric_means_features/enhancers_GMcounts.bed

### Merge HIV within Tad

