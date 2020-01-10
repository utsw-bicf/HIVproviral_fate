#! /bin/sh

#SBATCH --job-name=annotation
#SBATCH --partition=super
#SBATCH --output=annotation.%j.out
#SBATCH --error=annotation.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### Make bed file of genes
convert2bed --input=gtf --output=bed </project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode.bed

#######################################################################
###This isn't as easy of a problem as first thought
### I need to make a decision tree


##### First, is this intergenic or genic
### Make bed files
perl /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Bhive/txt2bed.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_integrations.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_integrations.bed

### remove duplicate integration sites between reps
sort -u -k 1,3 -k 5,5 /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_integrations.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_integrations_u.bed

module load bedops/2.4.14
convert2bed --input=gtf --output=bed </project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf | awk '$8=="gene" {print $0}' | awk '$1 ~ /chr/ {print $0}' | grep "protein_coding" >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_gene.bed

### Find integration sites that overlap a gene
module load bedtools/2.26.0
bedtools intersect -wa -wb -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_integrations_u.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_gene.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic.bed

bedtools intersect -v -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_integrations_u.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_gene.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic.bed

##### Next, treat genic and intergenic differently
### Intergenic: Find the nearest 2 TSS sites or Find the up/down stream TSS and make a histogram of distances vs direction.
convert2bed --input=gtf --output=bed </project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf | awk '$8=="exon" {print $0}' | awk '$1 ~ /chr/ {print $0}' | grep "exon_number 1;" | grep "protein_coding" >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_TSS.bed

sort -u -k 1,2 /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_TSS.bed | sort -u -k 1,1 -k 3,3 >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_uTSS.bed

bedtools closest -iu -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_uTSS.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic_downstream.bed

bedtools closest -id -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_uTSS.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic_upstream.bed

cat /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic_downstream.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic_upstream.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic_both.bed

awk '$5==$14 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic_both.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic_same.bed

awk '$5!=$14 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic_both.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic_opposite.bed

bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_uTSS.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic_closest.bed

##### If genic, does it overlap more than 1 gene
### Yes, it does overlap more than 1 gene
awk 'n=x[$1, $2, $3]{print n"\n"$0;} {x[$1, $2, $3]=$0;}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_2plusgene.bed

# Intron vs Exon
convert2bed --input=gtf --output=bed </project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf | awk '$8=="exon" {print $0}' | awk '$1 ~ /chr/ {print $0}' | grep "protein_coding" >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_exon.bed

awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $5, $6, $7, $8}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_2plusgene.bed | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_2plusgene_r.bed

bedtools intersect -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_2plusgene_r.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_exon.bed) | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_2plusgene_exon.bed

comm -23 <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_2plusgene_r.bed) <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/genic_2plusgene_exon.bed) > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_2plusgene_notexon.bed



### No, it doesn't overlap more than 1 gene
comm -23 <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic.bed) <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_2plusgene.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_1gene.bed

# are the annotations the same as the nearest TSS
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_1gene.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_uTSS.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_1gene_closest.bed

awk 'BEGIN {FS="\t"}; $12==$22 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_1gene_closest.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_1gene_same_as_TSS.bed

awk 'BEGIN {FS="\t"}; $12!=$22 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_1gene_closest.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_1gene_diff_as_TSS.bed

# Intron vs Exon
bedtools intersect -wa -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_1gene.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_exon.bed) | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_1gene_exon.bed

comm -23 <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_1gene.bed) <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_1gene_exon.bed) > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_1gene_notexon.bed


##########
### Now make pie chart from the 5+ types of annotations

# 1 Intergenic and same direction as nearest TSS; answer 172
awk '$5==$14 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic_closest.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic_same_4pie.bed
# 2 Intergenic, downstream of TSS, opposite direction
awk '{if ($14=="+" && $2>$11 && $5!=$14) {print $0;} else if($14=="-" && $2<$10 && $5!=$14) print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic_closest.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic_convergent_4pie.bed
#3 Intergenic, upstream of TSS, opposite direction
awk '{if ($14=="+" && $2<$10 && $5!=$14) {print $0;} else if($14=="-" && $2>$11 && $5!=$14) print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic_closest.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/intergenic_divergent_4pie.bed
#4 Genic, same *AS ANNOTATED* (which may not be nearest TSS)
awk '$5==$14 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_1gene.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_same_4pie.bed
#5 Genic, convergent *AS ANNOTATED* (which may not be nearest TSS)
awk '$5!=$14 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_1gene.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_convergent_4pie.bed
##New
# 6 2+ genes, both same and convergent
awk '{print $1, $2, $3, $5, $14}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_2plusgene.bed | uniq -u >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_2plusgene_both_4pie.bed
# 7 2+ genes, same
comm -23 <(awk '{print $1, $2, $3, $5, $14}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_2plusgene.bed | sort) <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_2plusgene_both_4pie.bed) | uniq | awk '$4==$5 {print $0}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_2plusgene_same_4pie.bed 
# 8 2+ genes, convergent
comm -23 <(awk '{print $1, $2, $3, $5, $14}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_2plusgene.bed | sort) <(sort genic_2plusgene_both_4pie.bed) | uniq | awk '$4!=$5 {print $0}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/pie_chart_integration/genic_2plusgene_convergent_4pie.bed 


### Pie chart in R
#setwd("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate")
#slices <- c(172, 86, 96, 793, 822, 36)
#lbls <- c("Intergenic-Same", "Intergenic-Convergent", "Intergenic-Divergent", "Intragenic-Same", "Intragenic-Convergent", "Intragenic-Both")
#pdf("HIV_integration_pie.pdf")
#pie(slices, labels=lbls, main="B-HIVE integration")
#dev.off()


