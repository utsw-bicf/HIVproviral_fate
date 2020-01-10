#! /bin/bash

#SBATCH --job-name=hm
#SBATCH --partition=super
#SBATCH --output=hm.%j.out
#SBATCH --error=hm.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### This script makes a heat map of HIV sorted expression by each of the data types.

H3K27ac='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K27ac_pooled.fc_signal.bw'
H3K36me3='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/GSM1603209_GFP_H3K36me3_pooled.fc_signal.bw'
H3K79me3='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/GSM1603215_GFP_H3K79me3_pooled.fc_signal.bw'
H3K9me3='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/GSM1603227_GFP_H3K9me3_pooled.fc_signal.bw'
H3K4me1='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/GSM1603225_GFP_H3K4me1_pooled.fc_signal.bw'
H3K4me3='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K4me3_pooled.fc_signal.bw'
RNAse='/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/RNAseq.bw'
MNase='/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/callPeaksMACS/MNase-seq_DORSO_Not-deposited_pooled.fc_signal.bw'
DNase='/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/callPeaksMACS/DNAse-seq_ENCODE_GSM736501_pooled.fc_signal.bw'
TTseq='/archive/shared/DOrso_BICF/WashU_browser_data/ttseq_merged.bw'
TTseqFor='/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_fwd.bw'
TTseqRev='/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_rev.bw'
H3K27me3='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/SRR647929_H3K27me3_pooled.fc_signal.bw'
damID='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/call_peaks/damIDseq_signal.bw'

HIV='/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/HIV_SortedByExp.bed'

### Make heatmaps
module load deeptools/2.5.0.1

for his in H3K27ac H3K36me3 H3K79me3 H3K9me3 H3K4me1 H3K4me3 RNAse MNase DNase H3K27me3 damID; do
  Ifile=$(eval echo "\$$his")
  computeMatrix reference-point -S ${Ifile} \
    -R ${HIV} \
    -a 2000 \
    -b 2000 \
    --referencePoint center \
    --sortRegions keep \
    -p max/2 \
    -o ${his}.mat.gz

  plotHeatmap -m ${his}.mat.gz \
    --sortRegions no \
    --refPointLabel "HIV IS" \
    --samplesLabel ${his} \
    --colorList 'white,blue' \
    -out ${his}_IS.pdf
done

### TTseq separately
  computeMatrix reference-point -S ${TTseq} ${TTseqFor} ${TTseqRev} \
    -R ${HIV} \
    -a 2000 \
    -b 2000 \
    --referencePoint center \
    --sortRegions keep \
    -p max/2 \
    -o TTseq.mat.gz

  plotHeatmap -m TTseq.mat.gz \
    --sortRegions no \
    --refPointLabel "HIV IS" \
    --samplesLabel TTseq TTseqForward TTseqReverse \
    --colorList 'white,blue' \
    -out TTseq_IS.pdf
done
