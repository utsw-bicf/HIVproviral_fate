#! /bin/bash

########## This script makes a heatmap of the H3K27ac and MNase-seq
########## 2 kb around HIV IS and matches the ML output from Jeon
module load deeptools/2.5.0.1

##### start with 1559 insertions and divide by
##### High: 455, Intermediate: 753, Low:351
grep -v "chrUn" /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/HIV_SortedByExp.bed | grep -v "0.440206297003631" | grep -v "0.832659573217461" | sort -gr -k 5,5 | head -n 455 >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_high.bed
grep -v "chrUn" /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/HIV_SortedByExp.bed | grep -v "0.440206297003631" | grep -v "0.832659573217461" | sort -gr -k 5,5 | tail -n 351 >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_low.bed
grep -v "chrUn" /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/HIV_SortedByExp.bed | grep -v "0.440206297003631" | grep -v "0.832659573217461" | sort -gr -k 5,5 | head -n 1208 | tail -n 753 >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_intermediate.bed

##### Make heat map of H3K27ac
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K27ac_pooled.fc_signal.bw \
    -R /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_high.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_intermediate.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_low.bed \
    -a 2000 \
    -b 2000 \
    --referencePoint center \
    --sortRegions keep \
    -p max/2 \
    -o /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/H3K27ac_ML.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/H3K27ac_ML.mat.gz \
    --refPointLabel "HIV IS" \
    --samplesLabel "H3K27ac" \
    --regionsLabel "High" "Intermediate" "Low" \
    --colors "#F8766D" "#00BA38" "#619CFF" \
    -out H3K27ac_IS_ML_profileplot.pdf


##### Make heat map of MNase-seq
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/callPeaksMACS/MNase-seq_DORSO_Not-deposited_pooled.fc_signal.bw \
    -R /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_high.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_intermediate.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_low.bed \
    -a 2000 \
    -b 2000 \
    --referencePoint center \
    --sortRegions keep \
    -p max/2 \
    -o /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/MNase_ML.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/MNase_ML.mat.gz \
    --refPointLabel "HIV IS" \
    --samplesLabel "MNase" \
    --regionsLabel "High" "Intermediate" "Low" \
    --colors "#F8766D" "#00BA38" "#619CFF" \
    -out MNase_IS_ML_profileplot.pdf


########################################################################
########## Do other marks
##### Make heat map of H3K4me3
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K4me3_pooled.fc_signal.bw \
    -R /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_high.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_intermediate.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_low.bed \
    -a 2000 \
    -b 2000 \
    --referencePoint center \
    --sortRegions keep \
    -p max/2 \
    -o /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/H3K4me3_ML.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/H3K4me3_ML.mat.gz \
    --refPointLabel "HIV IS" \
    --samplesLabel "H3K4me3" \
    --regionsLabel "High" "Intermediate" "Low" \
    --colors "#F8766D" "#00BA38" "#619CFF" \
    -out H3K4me3_IS_ML_profileplot.pdf


##### Make heat map of H3K4me1
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/GSM1603225_GFP_H3K4me1_pooled.fc_signal.bw \
    -R /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_high.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_intermediate.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_low.bed \
    -a 2000 \
    -b 2000 \
    --referencePoint center \
    --sortRegions keep \
    -p max/2 \
    -o /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/H3K4me1_ML.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/H3K4me1_ML.mat.gz \
    --refPointLabel "HIV IS" \
    --samplesLabel "H3K4me1" \
    --regionsLabel "High" "Intermediate" "Low" \
    --colors "#F8766D" "#00BA38" "#619CFF" \
    -out H3K4me1_IS_ML_profileplot.pdf


##### Make heat map of H3K36me3
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/GSM1603209_GFP_H3K36me3_pooled.fc_signal.bw \
    -R /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_high.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_intermediate.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_low.bed \
    -a 2000 \
    -b 2000 \
    --referencePoint center \
    --sortRegions keep \
    -p max/2 \
    -o /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/H3K36me3_ML.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/H3K36me3_ML.mat.gz \
    --refPointLabel "HIV IS" \
    --samplesLabel "H3K36me3" \
    --regionsLabel "High" "Intermediate" "Low" \
    --colors "#F8766D" "#00BA38" "#619CFF" \
    -out H3K36me3_IS_ML_profileplot.pdf


##### Make heat map of H3K79me3
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/GSM1603215_GFP_H3K79me3_pooled.fc_signal.bw \
    -R /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_high.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_intermediate.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_low.bed \
    -a 2000 \
    -b 2000 \
    --referencePoint center \
    --sortRegions keep \
    -p max/2 \
    -o /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/H3K79me3_ML.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/H3K79me3_ML.mat.gz \
    --refPointLabel "HIV IS" \
    --samplesLabel "H3K79me3" \
    --regionsLabel "High" "Intermediate" "Low" \
    --colors "#F8766D" "#00BA38" "#619CFF" \
    -out H3K79me3_IS_ML_profileplot.pdf


##### Make heat map of H3K9me3
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/GSM1603227_GFP_H3K9me3_pooled.fc_signal.bw \
    -R /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_high.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_intermediate.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_low.bed \
    -a 2000 \
    -b 2000 \
    --referencePoint center \
    --sortRegions keep \
    -p max/2 \
    -o /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/H3K9me3_ML.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/H3K9me3_ML.mat.gz \
    --refPointLabel "HIV IS" \
    --samplesLabel "H3K9me3" \
    --regionsLabel "High" "Intermediate" "Low" \
    --colors "#F8766D" "#00BA38" "#619CFF" \
    -out H3K9me3_IS_ML_profileplot.pdf



##### Make heat map of H3K9me3
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/SRR647929_H3K27me3_pooled.fc_signal.bw \
    -R /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_high.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_intermediate.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_low.bed \
    -a 2000 \
    -b 2000 \
    --referencePoint center \
    --sortRegions keep \
    -p max/2 \
    -o /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/H3K27me3_ML.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/H3K27me3_ML.mat.gz \
    --refPointLabel "HIV IS" \
    --samplesLabel "H3K27me3" \
    --regionsLabel "High" "Intermediate" "Low" \
    --colors "#F8766D" "#00BA38" "#619CFF" \
    -out H3K27me3_IS_ML_profileplot.pdf


##### Make heat map of RNAseq
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/RNAseq.bw \
    -R /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_high.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_intermediate.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_low.bed \
    -a 2000 \
    -b 2000 \
    --referencePoint center \
    --sortRegions keep \
    -p max/2 \
    -o /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/RNAseq_ML.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/RNAseq_ML.mat.gz \
    --refPointLabel "HIV IS" \
    --samplesLabel "RNAseq" \
    --regionsLabel "High" "Intermediate" "Low" \
    --colors "#F8766D" "#00BA38" "#619CFF" \
    -out RNAseq_IS_ML_profileplot.pdf


##### Make heat map of DNase
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/callPeaksMACS/DNAse-seq_ENCODE_GSM736501_pooled.fc_signal.bw \
    -R /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_high.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_intermediate.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_low.bed \
    -a 2000 \
    -b 2000 \
    --referencePoint center \
    --sortRegions keep \
    -p max/2 \
    -o /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/DNase_ML.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/DNase_ML.mat.gz \
    --refPointLabel "HIV IS" \
    --samplesLabel "DNase" \
    --regionsLabel "High" "Intermediate" "Low" \
    --colors "#F8766D" "#00BA38" "#619CFF" \
    -out DNase_IS_ML_profileplot.pdf


##### Make heat map of TTseq
computeMatrix reference-point -S /archive/shared/DOrso_BICF/WashU_browser_data/ttseq_merged.bw \
    -R /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_high.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_intermediate.bed /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_low.bed \
    -a 2000 \
    -b 2000 \
    --referencePoint center \
    --sortRegions keep \
    -p max/2 \
    -o /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/TTseq_ML.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/TTseq_ML.mat.gz \
    --refPointLabel "HIV IS" \
    --samplesLabel "TTseq" \
    --regionsLabel "High" "Intermediate" "Low" \
    --colors "#F8766D" "#00BA38" "#619CFF" \
    -out TTseq_IS_ML_profileplot.pdf
