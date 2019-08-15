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

convert2bed --input=gtf --output=bed </project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf | awk '$8=="exon" {print $0}' | awk '$1 ~ /chr/ {print $0}' | grep "protein_coding" >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_exon.bed

##### Turn HIV expression into a bed file
perl /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Bhive/exp_txt2bed.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression.bed

# Remove the duplicate lines between reps by selecting the highest number of reads
# Remove:chr4	79732674	79732675	GGCGCAAGTTCCTCCGACC	-	1	7	100	4969	2725	-0.260902490132872
# Remove: chr3	48249568	48249569	TGCACACCTAATCGT	+	1	6	100	2173	38	-1.75727612970371
# Leaves 1559 insertion/expression sites
grep -v "\-0.260902490132872" /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression.bed | grep -v "\-1.75727612970371" >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_u.bed

###############################################################################################
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles
########## Intergenic or Intragenic
bedtools intersect -wa -wb -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_u.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_gene.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic.bed

bedtools intersect -v -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/hiv_expression_u.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_gene.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intergenic.bed

###############################################################################################
########## Next, treat genic and intergenic differently
###############################################################################################

### Intergenic: Find the nearest TSS site, annotate it as such; Note: 1 cannot be annotated
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intergenic.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_uTSS.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intergenic_closestTSS.bed

### Intergenic: Find the nearest gene, annotate it as such
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intergenic.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_gene.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intergenic_closestGene.bed

###############################################################################################

### Genic: insertion in sites where there are 1 or 2+ genes
awk 'n=x[$1, $2, $3]{print n"\n"$0;} {x[$1, $2, $3]=$0;}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes.bed

comm -23 <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic.bed) <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_1gene.bed

###############################################################################################
########## Treat genic 1 gene and 2+ genes differently
###############################################################################################

########## Genic, 1 gene
### Genic, 1 gene: exon vs not exon
bedtools intersect -wa -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_1gene.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_exon.bed) | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_1gene_exon.bed

comm -23 <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_1gene.bed) <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_1gene_exon.bed) > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_1gene_notexon.bed

### Genic, 1 gene: Find the nearest TSS site, annotate it as such
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_1gene.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_uTSS.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_closestTSS.bed

### Genic, 1 gene: Find the nearest gene, annotate it as such
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_1gene.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_gene.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_closestGene.bed



########## Genic, 2+ gene
# Reduce the 2gene file of duplicates
awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $5, $6, $7, $8}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes.bed | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes_r.bed

### Genic, 2+ gene: exon vs not exon
bedtools intersect -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes_r.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_exon.bed) | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes_exon.bed

comm -23 <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes_r.bed) <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes_exon.bed) > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes_notexon.bed

### Genic, 2+ gene:: Find the nearest TSS site, annotate it as such
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes_r.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_uTSS.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes_closestTSS.bed

### Genic, 2+ gene:: Find the nearest gene, annotate it as such
bedtools closest -D a -a <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes_r.bed) -b <(bedtools sort -i /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/GRCh38_gencode_gene.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes_closestGene.bed



#####################################################################################################
### Find if closest TSS does not match closest Gene; 322 match and 34 don't, but 1 insertion is represented twice in Gene file
module load R/3.3.2-gccmkl
Rscript /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Bhive/hiv_expression_closestTSS_vs_closestGene.R


#####################################################################################################
### Now make pie chart from the 5+ types of annotations
rm /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt

# 1 Intergenic and same direction as nearest TSS
echo "# 1 Intergenic and same direction as nearest TSS" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt
awk '$5==$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intergenic_closestTSS.bed | wc -l >>/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt

# 2 Intergenic, downstream of TSS, opposite direction as nearest TSS
echo "# 2 Intergenic, downstream of TSS, opposite direction" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt
awk '{if ($17=="+" && $2>$14 && $5!=$17) {print $0;} else if($17=="-" && $2<$13 && $5!=$17) print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intergenic_closestTSS.bed | wc -l >>/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt

#3 Intergenic, upstream of TSS, opposite direction as nearest TSS
echo "#3 Intergenic, upstream of TSS, opposite direction" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt
awk '{if ($17=="+" && $2<$13 && $5!=$17) {print $0;} else if($17=="-" && $2>$14 && $5!=$17) print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intergenic_closestTSS.bed | wc -l >>/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt

#4 Genic, same
echo "#4 Genic, 1 gene" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt
awk '$5==$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_1gene.bed | wc -l >>/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt

#5 Genic, convergent
echo "#5 Genic, convergent" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt
awk '$5!=$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_1gene.bed | wc -l >>/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt

##New
# 6 2+ genes, both same and convergent
echo "# 6 2+ genes, both same and convergent" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt
awk '{print $1, $2, $3, $5, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes.bed | uniq -u | awk '{print $1, $2, $3, $4}' | sort | uniq | wc -l >>/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt
awk '{print $1, $2, $3, $5, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes.bed | uniq -u >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/delete.bed


# 7 2+ genes, same
echo "# 7 2+ genes, same" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt
comm -23 <(awk '{print $1, $2, $3, $5, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes.bed | sort) <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/delete.bed) | uniq | awk '$4==$5 {print $0}' | sort | uniq | wc -l >>/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt

# 8 2+ genes, convergent
echo "# 8 2+ genes, convergent" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt
comm -23 <(awk '{print $1, $2, $3, $5, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/hiv_expression_intragenic_2genes.bed | sort) <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/delete.bed) | uniq | awk '$4!=$5 {print $0}' | sort | uniq | wc -l >>/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/pie_counts.txt 

rm /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/delete.bed


### Pie chart in R
#setwd("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/")
#slices <- c(127, 71, 77, 607, 652, 24) #1,2,3,4+7,5+8,6
#lbls <- c("Intergenic-Same", "Intergenic-Convergent", "Intergenic-Divergent", "Intragenic-Same", "Intragenic-Convergent", "Intragenic-Overlapping")
#pdf("HIV_expression_pie_OLD.pdf")
#pie(slices, labels=lbls, main="B-HIVE integration")
#dev.off()

#setwd("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/")
#slices <- c(127, 71, 77, 581, 642, 60) #1,2,3,4,5,6+7+8
#lbls <- c("Intergenic-Same", "Intergenic-Convergent", "Intergenic-Divergent", "Intragenic-Same", "Intragenic-Convergent", "Intragenic-2genes")
#pdf("HIV_expression_pie.pdf")
#pie(slices, labels=lbls, main="B-HIVE integration")
#dev.off()



