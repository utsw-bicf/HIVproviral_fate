#! /bin/bash
### This graphic is all human genes vs H34me3, H3K27ac, and PolII heatmap

module load deeptools/2.5.0.1

### Human gene bed file
awk '{OFS="\t"}; $3=="gene" {print $1, $4, $5, $10, $6, $7}' /project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf | perl -p -e 's|"||g; s|;||g' >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/human_genes.bed

### Chip
H3K4me3='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K4me3_pooled.fc_signal.bw'
H3K27ac='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K27ac_pooled.fc_signal.bw'
PolII='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus_PolII/callPeaksMACS/PolII_pooled.fc_signal.bw'
DNAse='/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/callPeaksMACS/DNAse-seq_ENCODE_GSM736501_pooled.fc_signal.bw'
MNAseq='/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/callPeaksMACS/MNase-seq_DORSO_Not-deposited_pooled.pvalue_signal.bw'


### Make heatmaps, keep all and output sort order based on H3k4me3
computeMatrix reference-point -S ${H3K4me3} ${H3K27ac} ${PolII} ${DNAse} ${MNAseq} \
  -R human_genes.bed \
  -a 5000 \
  -b 5000 \
  --referencePoint TSS \
  --sortRegions descend \
  --sortUsing mean \
  --sortUsingSamples 1 \
  -p max/2 \
  --outFileSortedRegions humangenes_chip5_sort.bed \
  -o humangenes_chip5_sort.mat.gz

plotHeatmap -m humangenes_chip5_sort.mat.gz --refPointLabel "TSS" --samplesLabel H3K4me3 H3K27ac PolII DNAse MNAseq --colorList 'white,black' 'white,purple' 'white,green' 'white,blue' 'white,red' --sortRegions no -out humangenes_chip5.pdf

### Do it for RNAseq, keep the same order
awk '{OFS="\t"}; NR>1 {print $1,$2,$3,$4,$5,$6}' humangenes_chip5_sort.bed >humangenes_chip5_sort_reduced.bed
RNAseq='/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/RNAseq.bw'

computeMatrix reference-point -S ${RNAseq} \
  -R humangenes_chip5_sort_reduced.bed \
  -a 5000 \
  -b 5000 \
  --referencePoint TSS \
  --sortRegions keep \
  -p max/2 \
  -o humangenes_rnaseq_chipsort.mat.gz

plotHeatmap -m humangenes_rnaseq_chipsort.mat.gz --refPointLabel "TSS" --samplesLabel RNA --colorList 'white,red' --sortRegions no -out humangenes_rnaseq_chipsort.pdf
