#! /bin/bash
#SBATCH --job-name=plotprofile
#SBATCH --partition=super
#SBATCH --output=plotprofile.%j.out
#SBATCH --error=plotprofile.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

########## This script makes the plotprofile plots for figure 5B 
module load deeptools/2.5.0.1
module load bedtools/2.26.0

##### DNase-seq vs all enhancers and superenhancers
# enhancers
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/callPeaksMACS/DNAse-seq_ENCODE_GSM736501_pooled.fc_signal.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_Enhancer.bed \
  --referencePoint center \
  -a 2000 \
  -b 2000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/DNase_enhancers.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/DNase_enhancers.mat.gz \
  --samplesLabel DNase-seq \
  --refPointLabel center \
  --regionsLabel DNase-seq \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/DNase_enhancers.pdf

# super enhancers
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/callPeaksMACS/DNAse-seq_ENCODE_GSM736501_pooled.fc_signal.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed \
  --referencePoint center \
  -a 20000 \
  -b 20000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/DNase_superenhancers.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/DNase_superenhancers.mat.gz \
  --samplesLabel DNase-seq \
  --refPointLabel center \
  --regionsLabel DNase-seq \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/DNase_superenhancers.pdf



########################################################################
########################################################################
########################################################################
##### H3K27ac enhancers and superenhancers, divide by called with and call w/out
# enhancers
# To calculate w/ and w/out, will take all enhancers and remove overlap with particular mark
bedtools intersect -v -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_Enhancer.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/H3K27ac_enhancers.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/inputs/notH3K27ac_enhancers.bed

computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K27ac_pooled.fc_signal.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/H3K27ac_enhancers.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/inputs/notH3K27ac_enhancers.bed \
  --referencePoint center \
  -a 2000 \
  -b 2000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K27ac_enhancers.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K27ac_enhancers.mat.gz \
  --samplesLabel H3K27ac \
  --refPointLabel center \
  --regionsLabel "H3K27ac called peaks" "rest of peaks" \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K27ac_enhancers.pdf

# super enhancers
bedtools intersect -v -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_SEresults/H3K27ac_SE.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/inputs/notH3K27ac_superenhancers.bed

computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K27ac_pooled.fc_signal.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_SEresults/H3K27ac_SE.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/inputs/notH3K27ac_superenhancers.bed \
  --referencePoint center \
  -a 20000 \
  -b 20000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K27ac_superenhancers.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K27ac_superenhancers.mat.gz \
  --samplesLabel H3K27ac \
  --refPointLabel center \
  --regionsLabel "H3K27ac called peaks" "rest of peaks" \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K27ac_superenhancers.pdf


########################################################################
########################################################################
########################################################################
##### H3K4me1 enhancers and superenhancers, divide by called with and call w/out
# enhancers
# To calculate w/ and w/out, will take all enhancers and remove overlap with particular mark
bedtools intersect -v -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_Enhancer.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/H3K4me1_enhancers.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/inputs/notH3K4me1_enhancers.bed

computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/GSM1603225_GFP_H3K4me1_pooled.fc_signal.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/H3K4me1_enhancers.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/inputs/notH3K4me1_enhancers.bed \
  --referencePoint center \
  -a 2000 \
  -b 2000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K4me1_enhancers.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K4me1_enhancers.mat.gz \
  --samplesLabel H3K4me1 \
  --refPointLabel center \
  --regionsLabel "H3K4me1 called peaks" "rest of peaks" \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K4me1_enhancers.pdf

# super enhancers
bedtools intersect -v -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_SEresults/H3K4me1_SE.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/inputs/notH3K4me1_superenhancers.bed

computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/GSM1603225_GFP_H3K4me1_pooled.fc_signal.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_SEresults/H3K4me1_SE.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/inputs/notH3K4me1_superenhancers.bed \
  --referencePoint center \
  -a 20000 \
  -b 20000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K4me1_superenhancers.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K4me1_superenhancers.mat.gz \
  --samplesLabel H3K4me1 \
  --refPointLabel center \
  --regionsLabel "H3K4me1 called peaks" "rest of peaks" \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K4me1_superenhancers.pdf


########################################################################
########################################################################
########################################################################
##### H3K4me3 enhancers and superenhancers, divide by called with and call w/out
# enhancers
# To calculate w/ and w/out, will take all enhancers and remove overlap with particular mark
bedtools intersect -v -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_Enhancer.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/H3K4me3_enhancers.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/inputs/notH3K4me3_enhancers.bed

computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K4me3_pooled.fc_signal.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/H3K4me3_enhancers.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/inputs/notH3K4me3_enhancers.bed \
  --referencePoint center \
  -a 2000 \
  -b 2000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K4me3_enhancers.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K4me3_enhancers.mat.gz \
  --samplesLabel H3K4me3 \
  --refPointLabel center \
  --regionsLabel "H3K4me3 called peaks" "rest of peaks" \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K4me3_enhancers.pdf

# super enhancers
bedtools intersect -v -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_SEresults/H3K4me3_SE.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/inputs/notH3K4me3_superenhancers.bed

computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K4me3_pooled.fc_signal.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_SEresults/H3K4me3_SE.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/inputs/notH3K4me3_superenhancers.bed \
  --referencePoint center \
  -a 20000 \
  -b 20000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K4me3_superenhancers.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K4me3_superenhancers.mat.gz \
  --samplesLabel H3K4me3 \
  --refPointLabel center \
  --regionsLabel "H3K4me3 called peaks" "rest of peaks" \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/H3K4me3_superenhancers.pdf



##### TT-seq vs all enhancers and superenhancers
# enhancers
computeMatrix reference-point -S /archive/shared/DOrso_BICF/WashU_browser_data/ttseq_merged.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_Enhancer.bed \
  --referencePoint center \
  -a 2000 \
  -b 2000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/TTseq_enhancers.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/TTseq_enhancers.mat.gz \
  --samplesLabel TT-seq \
  --refPointLabel center \
  --regionsLabel TT-seq \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/TTseq_enhancers.pdf

# super enhancers
computeMatrix reference-point -S /archive/shared/DOrso_BICF/WashU_browser_data/ttseq_merged.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed \
  --referencePoint center \
  -a 20000 \
  -b 20000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/TTseq_superenhancers.mat.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/TTseq_superenhancers.mat.gz \
  --samplesLabel TT-seq \
  --refPointLabel center \
  --regionsLabel TT-seq \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/profile_fig5B/TTseq_superenhancers.pdf

