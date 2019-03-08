#! /bin/bash
### This script annotates the intertions of the BHIVE expression data
### only 1572 of the 2010 sites are consistant between the integration data and expression data

### Load programs
module load bedops/2.4.14
module load bedtools/2.26.0

########## Make files into bed files
### Make bed file from gtf file
convert2bed --input=gtf --output=bed </project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode.bed

convert2bed --input=gtf --output=bed </project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf | awk '$8=="gene" {print $0}' | awk '$1 ~ /chr/ {print $0}' | grep "protein_coding" >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_gene.bed

convert2bed --input=gtf --output=bed </project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf | awk '$8=="exon" {print $0}' | awk '$1 ~ /chr/ {print $0}' | grep "exon_number 1;" | grep "protein_coding" >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_TSS.bed
sort -u -k 1,2 /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_TSS.bed | sort -u -k 1,1 -k 3,3 >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_uTSS.bed

##### Turn HIV expression into a bed file
perl /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Bhive/exp_txt2bed.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression.bed

# Remove the duplicate lines between reps by selecting the highest mapping quality
# Remove:chr3	48249568	48249569	TGCACACCTAATCGT	+	2	35956	26	106008	181932	0.0794704891172975
# Remove:chr4	79732674	79732675	GGCGCAAGTTCCTCCGACC	-	2	39441	15	1056969	3840998	0.405281867316955
# Leaves 1570 insertion/expression sites
grep -v "0.0794704891172975" /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression.bed | grep -v "0.405281867316955" >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_u.bed

###############################################################################################

########## Intergenic or Intragenic
bedtools intersect -wa -wb -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_u.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_gene.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic.bed

bedtools intersect -v -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_integrations_u.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_gene.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intergenic.bed

###############################################################################################
########## Next, treat genic and intergenic differently
###############################################################################################

### Intergenic: Find the nearest TSS site, annotate it as such
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intergenic.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_uTSS.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intergenic_closestTSS.bed

### Intergenic: Find the nearest gene, annotate it as such
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intergenic.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_gene.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intergenic_closestGene.bed

### Find if closest TSS does not match closest Gene; 322 match and 34 don't, but 1 insertion is represented twice in Gene file
module load R/3.3.2-gccmkl
Rscript /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Bhive/hiv_expression_closestTSS_vs_closestGene.R
