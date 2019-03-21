#! /bin/bash
### This script makes a figure that is a heatmap of
### expression insertions by location (from pie chart)
### by chipseq analysis: H3K4me3, H3K27ac, PolII, chromatin accessibility
### Using deeptools

########## Make bed files from pie chart
# Make a delete.bed
awk '{print $1, $2, $3, $4, $11, $5, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic_2genes.bed | uniq -u >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/delete.bed

# 1
awk '$5==$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intergenic_closestTSS.bed | awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $11, $5}' | sort -nr -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number1.bed

# 2
awk '{if ($17=="+" && $2>$14 && $5!=$17) {print $0;} else if($17=="-" && $2<$13 && $5!=$17) print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intergenic_closestTSS.bed | awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $11, $5}' | sort -nr -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number2.bed

# 3
awk '{if ($17=="+" && $2<$13 && $5!=$17) {print $0;} else if($17=="-" && $2>$14 && $5!=$17) print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intergenic_closestTSS.bed | awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $11, $5}' | sort -nr -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number3.bed

# 4
awk '$5==$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic_1gene.bed | awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $11, $5}' | sort -nr -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number4.bed

# 5
awk '$5!=$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic_1gene.bed | awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $11, $5}' | sort -nr -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number5.bed

# 6
awk '{print $1, $2, $3, $4, $11, $5, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic_2genes.bed | uniq -u | awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $11, $5}' | sort | uniq | sort -nr -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number6.bed

#7
comm -23 <(awk '{print $1, $2, $3, $4, $11, $5, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic_2genes.bed | sort) <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/delete.bed) | uniq | awk 'BEGIN {OFS="\t"}; $6==$7 {print $1, $2, $3, $4, $5, $6}' | sort | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number7.bed

#8
comm -23 <(awk '{print $1, $2, $3, $4, $11, $5, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic_2genes.bed | sort) <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/delete.bed) | uniq | awk 'BEGIN {OFS="\t"}; $6!=$7 {print $1, $2, $3, $4, $5, $6}' | sort | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number8.bed 

# 4+7
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number4.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number7.bed | sort -nr -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number4_7.bed

# 5+8
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number5.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number8.bed | sort -nr -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number5_8.bed

# Remove delete.bed
rm /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/delete.bed




########## Location of bw files
H3K4me3='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K4me3_pooled.fc_signal.bw'
H3K27ac='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K27ac_pooled.fc_signal.bw'
PolII='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus_PolII/callPeaksMACS/PolII_pooled.fc_signal.bw'
#DNAse='/project/BICF/BICF_Core/s185797/testing_ATACseq4DNAseseq/atacseq_analysis/workflow/output/callPeaksMACS/DNAse-seq_ENCODE_GSM736501_pooled.fc_signal.bw'
DNAse='/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/callPeaksMACS/DNAse-seq_ENCODE_GSM736501_pooled.fc_signal.bw'

########## Heatmaps
module load deeptools/2.5.0.1
computeMatrix reference-point -S ${H3K4me3} ${H3K27ac} ${PolII} ${DNAse} \
  -R number4_7.bed number5_8.bed number1.bed number2.bed number3.bed \
  -a 5000 \
  -b 5000 \
  --referencePoint center \
  --sortRegions no \
  --skipZeros \
  -p max/2 \
  -o 1D.mat.gz

plotHeatmap -m 1D.mat.gz --sortRegions no --refPointLabel "HIV IS" --samplesLabel H3K4me3 H3K27ac PolII DNase --regionsLabel Same Opposite Same C D --colorList 'white,black' 'white,purple' 'white,green' 'white,blue' -out 1D_IS.pdf


####################################################################################################
### There was no pattern with the above heatmaps 
### So I will create a heatmap of the nearest TSS genes
####################################################################################################
# delete
awk '{print $1, $2, $3, $4, $11, $5, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic_2genes.bed | uniq -u >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/delete.bed

# 1
awk '$5==$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intergenic_closestTSS.bed | awk 'BEGIN {OFS="\t"}; {print $12, $13, $14, $15, $11, $17}' | sort -nr -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number1_genes.bed

# 2
awk '{if ($17=="+" && $2>$14 && $5!=$17) {print $0;} else if($17=="-" && $2<$13 && $5!=$17) print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intergenic_closestTSS.bed | awk 'BEGIN {OFS="\t"}; {print $12, $13, $14, $15, $11, $17}' | sort -nr -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number2_genes.bed

# 3
awk '{if ($17=="+" && $2<$13 && $5!=$17) {print $0;} else if($17=="-" && $2>$14 && $5!=$17) print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intergenic_closestTSS.bed | awk 'BEGIN {OFS="\t"}; {print $12, $13, $14, $15, $11, $17}' | sort -nr -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number3_genes.bed

# 4
awk '$5==$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic_1gene.bed | awk 'BEGIN {OFS="\t"}; {print $12, $13, $14, $15, $11, $17}' | sort -nr -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number4_genes.bed

# 5
awk '$5!=$17 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic_1gene.bed | awk 'BEGIN {OFS="\t"}; {print $12, $13, $14, $15, $11, $17}' | sort -nr -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number5_genes.bed

# 6
awk '{print $1, $2, $3, $4, $11, $5, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic_2genes.bed | uniq -u | awk '{print $5}' | sort | uniq | sort -nr -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number6_delete.bed
grep -f /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number6_delete.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic_2genes.bed | awk 'BEGIN {OFS="\t"}; {if ($5==$17) print $12, $13, $14, $15, $11, $17;}' | uniq | sort -nr -k 5,5 >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number6_genes_same.bed
grep -f /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number6_delete.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic_2genes.bed | awk 'BEGIN {OFS="\t"}; {if ($5!=$17) print $12, $13, $14, $15, $11, $17;}' | uniq | sort -nr -k 5,5 >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number6_genes_opposite.bed

#7
comm -23 <(awk '{print $1, $2, $3, $4, $11, $5, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic_2genes.bed | sort) <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/delete.bed) | uniq | awk 'BEGIN {OFS="\t"}; $6==$7 {print $5}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number7_delete.bed
grep -f /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number7_delete.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic_2genes.bed | awk 'BEGIN {OFS="\t"}; {print $12, $13, $14, $15, $11, $17}' | uniq | sort -nr -k 5,5 >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number7_genes.bed

#8
comm -23 <(awk '{print $1, $2, $3, $4, $11, $5, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic_2genes.bed | sort) <(sort /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/delete.bed) | uniq | awk 'BEGIN {OFS="\t"}; $6!=$7 {print $5}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number8_delete.bed
grep -f /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number8_delete.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated/hiv_expression_intragenic_2genes.bed | awk 'BEGIN {OFS="\t"}; {print $12, $13, $14, $15, $11, $17}' | uniq | sort -nr -k 5,5 >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number8_genes.bed

rm /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/*delete*


########## Location of bw files
H3K4me3='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K4me3_pooled.fc_signal.bw'
H3K27ac='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K27ac_pooled.fc_signal.bw'
PolII='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus_PolII/callPeaksMACS/PolII_pooled.fc_signal.bw'
#DNAse='/project/BICF/BICF_Core/s185797/testing_ATACseq4DNAseseq/atacseq_analysis/workflow/output/callPeaksMACS/DNAse-seq_ENCODE_GSM736501_pooled.fc_signal.bw'
DNAse='/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/callPeaksMACS/DNAse-seq_ENCODE_GSM736501_pooled.fc_signal.bw'

########## Heatmaps
module load deeptools/2.5.0.1
computeMatrix scale-regions -S ${H3K4me3} ${H3K27ac} ${PolII} ${DNAse} \
  -R number4_genes.bed number5_genes.bed number1_genes.bed number2_genes.bed number3_genes.bed \
  -a 5000 \
  -b 5000 \
  -m 5000 \
  --sortRegions no \
  --skipZeros \
  -p max/2 \
  -o 1D_genes.mat.gz

plotHeatmap -m 1D_genes.mat.gz --sortRegions no --samplesLabel H3K4me3 H3K27ac PolII DNase --regionsLabel Same Opposite Same C D --colorList 'white,black' 'white,purple' 'white,green' 'white,blue' -out 1D_genes.pdf


####################################################################################################
### Still nothing; add RNAseq to the figures
####################################################################################################
### RNAseq bam to bw
module load samtools 
samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/RNAseq.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq_reg/workflow/output_PE/Emily_jurkat_*.dedup.bam
samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/RNAseq.bam
bamCoverage -p max/2 -b /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/RNAseq.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/RNAseq.bw

########## Location of bw files
H3K4me3='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K4me3_pooled.fc_signal.bw'
H3K27ac='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K27ac_pooled.fc_signal.bw'
PolII='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus_PolII/callPeaksMACS/PolII_pooled.fc_signal.bw'
#DNAse='/project/BICF/BICF_Core/s185797/testing_ATACseq4DNAseseq/atacseq_analysis/workflow/output/callPeaksMACS/DNAse-seq_ENCODE_GSM736501_pooled.fc_signal.bw'
DNAse='/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/callPeaksMACS/DNAse-seq_ENCODE_GSM736501_pooled.fc_signal.bw'
RNAseq='/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/RNAseq.bw'

########## Heatmaps
computeMatrix reference-point -S ${H3K4me3} ${H3K27ac} ${PolII} ${DNAse} ${RNAseq} \
  -R number4_7.bed number5_8.bed number1.bed number2.bed number3.bed \
  -a 5000 \
  -b 5000 \
  --referencePoint center \
  --sortRegions no \
  --skipZeros \
  -p max/2 \
  -o 1D_ChipRNA.mat.gz

plotHeatmap -m 1D_ChipRNA.mat.gz --sortRegions no --refPointLabel "HIV IS" --samplesLabel H3K4me3 H3K27ac PolII DNase RNA --regionsLabel Same Opposite Same C D --colorList 'white,black' 'white,purple' 'white,green' 'white,blue' 'white,red' -out 1D_ChipRNA_IS.pdf

computeMatrix scale-regions -S ${H3K4me3} ${H3K27ac} ${PolII} ${DNAse} ${RNAseq} \
  -R number4_genes.bed number5_genes.bed number1_genes.bed number2_genes.bed number3_genes.bed \
  -a 5000 \
  -b 5000 \
  -m 5000 \
  --sortRegions no \
  --skipZeros \
  -p max/2 \
  -o 1D_ChipRNA_genes.mat.gz

plotHeatmap -m 1D_ChipRNA_genes.mat.gz --sortRegions no --samplesLabel H3K4me3 H3K27ac PolII DNase RNA --regionsLabel Same Opposite Same C D --colorList 'white,black' 'white,purple' 'white,green' 'white,blue' 'white,red' -out 1D_ChipRNA_genes.pdf

### Heatmaps look bad because the scale of chipseq vs rnaseq is different; so just run rnaseq
computeMatrix reference-point -S ${RNAseq} \
  -R number4_7.bed number5_8.bed number1.bed number2.bed number3.bed \
  -a 5000 \
  -b 5000 \
  --referencePoint center \
  --sortRegions no \
  --skipZeros \
  -p max/2 \
  -o 1D_RNA.mat.gz

plotHeatmap -m 1D_RNA.mat.gz --sortRegions no --refPointLabel "HIV IS" --samplesLabel RNA --regionsLabel Same Opposite Same C D --colorList 'white,red' -out 1D_RNA_IS.pdf

computeMatrix scale-regions -S ${RNAseq} \
  -R number4_genes.bed number5_genes.bed number1_genes.bed number2_genes.bed number3_genes.bed \
  -a 5000 \
  -b 5000 \
  -m 5000 \
  --sortRegions no \
  --skipZeros \
  -p max/2 \
  -o 1D_RNA_genes.mat.gz

plotHeatmap -m 1D_RNA_genes.mat.gz --sortRegions no --samplesLabel RNA --regionsLabel Same Opposite Same C D --colorList 'white,red' -out 1D_RNA_genes.pdf

### Add mnase seq also
H3K4me3='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K4me3_pooled.fc_signal.bw'
H3K27ac='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K27ac_pooled.fc_signal.bw'
PolII='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus_PolII/callPeaksMACS/PolII_pooled.fc_signal.bw'
#DNAse='/project/BICF/BICF_Core/s185797/testing_ATACseq4DNAseseq/atacseq_analysis/workflow/output/callPeaksMACS/DNAse-seq_ENCODE_GSM736501_pooled.fc_signal.bw'
DNAse='/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/callPeaksMACS/DNAse-seq_ENCODE_GSM736501_pooled.fc_signal.bw'
MNAseq='/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/callPeaksMACS/MNase-seq_DORSO_Not-deposited_pooled.pvalue_signal.bw'

########## Heatmaps
computeMatrix reference-point -S ${H3K4me3} ${H3K27ac} ${PolII} ${DNAse} ${MNAseq} \
  -R number4_7.bed number5_8.bed number1.bed number2.bed number3.bed \
  -a 5000 \
  -b 5000 \
  --referencePoint center \
  --sortRegions no \
  --skipZeros \
  -p max/2 \
  -o 1D_Chip5.mat.gz

plotHeatmap -m 1D_Chip5.mat.gz --sortRegions no --refPointLabel "HIV IS" --samplesLabel H3K4me3 H3K27ac PolII DNase MNase --regionsLabel Same Opposite Same C D --colorList 'white,black' 'white,purple' 'white,green' 'white,blue' 'white,red' -out 1D_Chip5_IS.pdf

computeMatrix scale-regions -S ${H3K4me3} ${H3K27ac} ${PolII} ${DNAse} ${MNAseq} \
  -R number4_genes.bed number5_genes.bed number1_genes.bed number2_genes.bed number3_genes.bed \
  -a 5000 \
  -b 5000 \
  -m 5000 \
  --sortRegions no \
  --skipZeros \
  -p max/2 \
  -o 1D_Chip5_genes.mat.gz

plotHeatmap -m 1D_Chip5_genes.mat.gz --sortRegions no --samplesLabel H3K4me3 H3K27ac PolII DNase MNase --regionsLabel Same Opposite Same C D --colorList 'white,black' 'white,purple' 'white,green' 'white,blue' 'white,red' -out 1D_Chip5_genes.pdf

