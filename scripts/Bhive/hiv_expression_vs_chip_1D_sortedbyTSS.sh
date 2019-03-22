#! /bin/bash
#### This script makes a heatmap of each HIV expression cluster vs chip and RNAseq, sorted by TSS distance

module load deeptools/2.5.0.1

### Make HIV expression bed files
# 1
awk -F "\t" 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $23, $5}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed | sort -n -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number1_disTSS.bed

# 2
awk -F "\t" 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $23, $5}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed | sort -n -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number2_disTSS.bed

# 3
awk -F "\t" 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $23, $5}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed | sort -n -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number3_disTSS.bed

# 4
awk -F "\t" 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $23, $5}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed | sort -n -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number4_disTSS.bed

# 5
awk -F "\t" 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $23, $5}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed | sort -n -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number5_disTSS.bed



### CHIPs
H3K4me3='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K4me3_pooled.fc_signal.bw'
H3K27ac='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K27ac_pooled.fc_signal.bw'
PolII='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus_PolII/callPeaksMACS/PolII_pooled.fc_signal.bw'
DNAse='/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/callPeaksMACS/DNAse-seq_ENCODE_GSM736501_pooled.fc_signal.bw'
MNAseq='/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/callPeaksMACS/MNase-seq_DORSO_Not-deposited_pooled.pvalue_signal.bw'


########## Heatmaps
computeMatrix reference-point -S ${H3K4me3} ${H3K27ac} ${PolII} ${DNAse} ${MNAseq} \
  -R number4_disTSS.bed number5_disTSS.bed number1_disTSS.bed number2_disTSS.bed number3_disTSS.bed \
  -a 5000 \
  -b 5000 \
  --referencePoint center \
  --sortRegions keep \
  -p max/2 \
  --outFileSortedRegions 1D_Chip5_sortbydisTSS.bed \
  -o 1D_Chip5_sortbydisTSS.mat.gz

plotHeatmap -m 1D_Chip5_sortbydisTSS.mat.gz --sortRegions no --refPointLabel "HIV IS" --samplesLabel H3K4me3 H3K27ac PolII DNase MNase --regionsLabel Same Opposite Same C D --colorList 'white,black' 'white,purple' 'white,green' 'white,blue' 'white,red' -out 1D_Chip5_IS_sortbydisTSS.pdf


########## By genes
### Make bed files of HIV IS in genes
# 1
awk -F "\t" 'BEGIN {OFS="\t"}; {print $12, $13, $14, $4, $23, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed | sort -n -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number1_genes_disTSS.bed

# 2
awk -F "\t" 'BEGIN {OFS="\t"}; {print $12, $13, $14, $4, $23, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed | sort -n -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number2_genes_disTSS.bed

# 3
awk -F "\t" 'BEGIN {OFS="\t"}; {print $12, $13, $14, $4, $23, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed | sort -n -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number3_genes_disTSS.bed

# 4
awk -F "\t" 'BEGIN {OFS="\t"}; {print $12, $13, $14, $4, $23, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed | sort -n -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number4_genes_disTSS.bed

# 5
awk -F "\t" 'BEGIN {OFS="\t"}; {print $12, $13, $14, $4, $23, $17}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed | sort -n -k 5,5 > /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/number5_genes_disTSS.bed

computeMatrix scale-regions -S ${H3K4me3} ${H3K27ac} ${PolII} ${DNAse} ${MNAseq} \
  -R number4_genes_disTSS.bed number5_genes_disTSS.bed number1_genes_disTSS.bed number2_genes_disTSS.bed number3_genes_disTSS.bed \
  -a 5000 \
  -b 5000 \
  -m 5000 \
  --sortRegions keep \
  -p max/2 \
  --outFileSortedRegions 1D_Chip5_genes_sortbydisTSS.bed \
  -o 1D_Chip5_genes_sortbydisTSS.mat.gz

plotHeatmap -m 1D_Chip5_genes_sortbydisTSS.mat.gz --sortRegions no --samplesLabel H3K4me3 H3K27ac PolII DNase MNase --regionsLabel Same Opposite Same C D --colorList 'white,black' 'white,purple' 'white,green' 'white,blue' 'white,red' -out 1D_Chip5_genes_sortbydisTSS.pdf


########## RNA
RNAseq='/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/RNAseq.bw'
computeMatrix reference-point -S ${RNAseq} \
  -R number4_disTSS.bed number5_disTSS.bed number1_disTSS.bed number2_disTSS.bed number3_disTSS.bed \
  -a 5000 \
  -b 5000 \
  --referencePoint center \
  --sortRegions keep \
  -p max/2 \
  --outFileSortedRegions 1D_RNA_sortbydisTSS.bed \
  -o 1D_RNA_sortbydisTSS.mat.gz

plotHeatmap -m 1D_RNA_sortbydisTSS.mat.gz --sortRegions no --refPointLabel "HIV IS" --samplesLabel RNA --regionsLabel Same Opposite Same C D --colorList 'white,red' -out 1D_RNA_IS_sortbydisTSS.pdf

computeMatrix scale-regions -S ${RNAseq} \
  -R number4_genes_disTSS.bed number5_genes_disTSS.bed number1_genes_disTSS.bed number2_genes_disTSS.bed number3_genes_disTSS.bed \
  -a 5000 \
  -b 5000 \
  -m 5000 \
  --sortRegions keep \
  -p max/2 \
  --outFileSortedRegions 1D_RNA_genes_sortbydisTSS.bed \
  -o 1D_RNA_genes_sortbydisTSS.mat.gz

plotHeatmap -m 1D_RNA_genes_sortbydisTSS.mat.gz --sortRegions no --samplesLabel RNA --regionsLabel Same Opposite Same C D --colorList 'white,red' -out 1D_RNA_genes_sortbydisTSS.pdf
