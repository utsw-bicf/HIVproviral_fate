#!/bin/bash
### This script annotates the distance to the nearest upstream TSS site.
### First I will do it Ivan's way: x=log10(distance to TSS) y=(HIV expression) but this doesn't take into account whether it's in a gene.
### Then I'll do it my way. Look at HIV expression of genic insertions and relate it to RNAseq expression
### Then take the intergenic HIV expression and see if there are enhancers upstream. An enhancer will be defined as TT-seq that overlaps with H3K27Ac that isn't within 2000 bp of TSS. Then compare TT-seq expression with HIV expression

########## Ivan's way: x=log10(distance to TSS) y=(HIV expression)
##### Turn HIV expression into a bed file
perl /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Bhive/exp_txt2bed.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression.bed

### don't remove duplicate integration sites between reps, since there are only 2 and expression doesn't appear uniform
##### Annotate the up and downstream gene
### Note: get uTSS file from bhive_annotate.sh
module load bedtools/2.26.0

bedtools closest -iu -D b -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_uTSS.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_proximalTSS.bed

### Divide by orientation: same or opposite
awk '$5==$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_proximalTSS.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_proximalTSS_same.bed
awk '$5!=$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_proximalTSS.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_proximalTSS_opposite.bed

### add newly calculated distance to bed files
awk 'BEGIN {OFS="\t"} {if ($17=="+") {print $0, $2-$13;} else if($17=="-") print $0, $14-$3}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_proximalTSS_same.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_proximalTSS_same_plus_dis.bed

awk 'BEGIN {OFS="\t"} {if ($17=="+") {print $0, $2-$13;} else if($17=="-") print $0, $14-$3}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_proximalTSS_opposite.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_proximalTSS_opposite_plus_dis.bed

##############################################################################################################
########## My way: HIV expression of genic insertions and relate it to RNAseq expression
### Find only that overlaps with a gene
### Note : get GRCh38_gencode_gene.bed from bhive_annotate.sh
bedtools intersect -wa -wb -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_gene.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_genic.bed

### Divide by orientation: same or opposite
awk '$5==$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_genic.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_genic_same.bed
awk '$5!=$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_genic.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_genic_opposite.bed

### add newly calculated distance to bed files
awk 'BEGIN {OFS="\t"} {if ($17=="+") {print $0, $2-$13;} else if($17=="-") print $0, $14-$3}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_genic_same.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_genic_same_plus_dis.bed

awk 'BEGIN {OFS="\t"} {if ($17=="+") {print $0, $2-$13;} else if($17=="-") print $0, $14-$3}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_genic_opposite.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_genic_opposite_plus_dis.bed

##############################################################################################################
########## My way: HIV expression of exonic insertions and relate it to RNAseq expression
### Note : get GRCh38_gencode_exon.bed from bhive_annotate.sh
bedtools intersect -wa -wb -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_exon.bed) | sort -u -k 1,5 >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_exon.bed

### Divide by orientation: same or opposite
awk '$5==$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_exon.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_exon_same.bed
awk '$5!=$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_exon.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_exon_opposite.bed

### add newly calculated distance to bed files
awk 'BEGIN {OFS="\t"} {if ($17=="+") {print $0, $2-$13;} else if($17=="-") print $0, $14-$3}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_exon_same.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_exon_same_plus_dis.bed

awk 'BEGIN {OFS="\t"} {if ($17=="+") {print $0, $2-$13;} else if($17=="-") print $0, $14-$3}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_exon_opposite.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_exon_opposite_plus_dis.bed

##############################################################################################################
### Make scatter plot in R
module load R/3.3.2-gccmkl
Rscript /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Bhive/BHIVE_expression_scatterplots.R
