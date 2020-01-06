#! /bin/bash
### This plots the Pearson Correlation for chipseq data
### Makes plot pdf
### Make pca

module load deeptools/2.5.0.1

###Make Pearson Correlation
plotCorrelation \
  -in /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/work/bf/a248cb17864b87cd69c398cbd9bf7b/sample_mbs.npz \
  --corMethod pearson \
  --skipZero \
  --plotTitle "Pearson Correlation of Read Counts" \
  --whatToPlot heatmap \
  --colorMap RdYlBu \
  --plotNumbers \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/heatmap_PearsonCorr.pdf

###Make Spearman Correlation:pdf
plotCorrelation \
  -in /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/work/bf/a248cb17864b87cd69c398cbd9bf7b/sample_mbs.npz \
  --corMethod spearman \
  --skipZero \
  --plotTitle "Spearman Correlation of Read Counts" \
  --whatToPlot heatmap \
  --colorMap RdYlBu \
  --plotNumbers \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/heatmap_SpearmanCorr.pdf

###Make PCA of readcounts
plotPCA \
  -in /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/work/bf/a248cb17864b87cd69c398cbd9bf7b/sample_mbs.npz \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/PCA_readCounts.pdf \
  -T "PCA of read counts"

###Make Pearson Correlation, removing outliers
plotCorrelation \
  -in /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/work/bf/a248cb17864b87cd69c398cbd9bf7b/sample_mbs.npz \
  --corMethod pearson \
  --skipZero \
  --plotTitle "Pearson Correlation of Read Counts" \
  --whatToPlot heatmap \
  --colorMap RdYlBu \
  --plotNumbers \
  --removeOutliers \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/heatmap_PearsonCorr_nooutliers.pdf

### Make individual pca plots
###Inputs
bam_loc='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/'
multiBamSummary bins \
  --bamfiles ${bam_loc}SRR2043612.filt.nodup.bam ${bam_loc}SRR1603648.filt.nodup.bam ${bam_loc}SRR1603649.filt.nodup.bam ${bam_loc}SRR1603652.filt.nodup.bam ${bam_loc}SRR1522118.filt.nodup.bam ${bam_loc}SRR061744.filt.nodup.bam ${bam_loc}SRR577484.filt.nodup.bam ${bam_loc}SRR074203.filt.nodup.bam ${bam_loc}SRR3722570.filt.nodup.bam ${bam_loc}SRR3722573.filt.nodup.bam ${bam_loc}SRR6010202.filt.nodup.bam ${bam_loc}SRR1057276.filt.nodup.bam ${bam_loc}GSM1603229.filt.nodup.bam \
  --extendReads 100 \
  --labels SRR2043612_input SRR1603648_input SRR1603649_input SRR1603652_input SRR1522118_input SRR061744_wce SRR577484_input SRR074203_input SRR3722570_WCE SRR3722573_DMSO_WCE SRR6010202_input SRR1057276_WCE GSM1603229_GFP_input \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/inputs_mbs.npz 

plotPCA \
  -in /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/inputs_mbs.npz \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/PCA_readCounts_inputs.pdf \
  -T "PCA of read counts - Inputs"

### Experiments
bam_loc='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/'
multiBamSummary bins \
  --bamfiles ${bam_loc}SRR2043614.filt.nodup.bam ${bam_loc}SRR1603646.filt.nodup.bam ${bam_loc}SRR1603650.filt.nodup.bam ${bam_loc}SRR1603654.filt.nodup.bam ${bam_loc}SRR1522115.filt.nodup.bam ${bam_loc}SRR061743.filt.nodup.bam ${bam_loc}SRR061747.filt.nodup.bam ${bam_loc}SRR647929.filt.nodup.bam ${bam_loc}SRR577482.filt.nodup.bam ${bam_loc}SRR577483.filt.nodup.bam ${bam_loc}SRR074201.filt.nodup.bam ${bam_loc}SRR074195.filt.nodup.bam ${bam_loc}SRR074196.filt.nodup.bam ${bam_loc}SRR3722566.filt.nodup.bam ${bam_loc}SRR3722577.filt.nodup.bam ${bam_loc}SRR6010201.filt.nodup.bam ${bam_loc}SRR2157609.filt.nodup.bam ${bam_loc}GSM1603225.filt.nodup.bam ${bam_loc}GSM1603213.filt.nodup.bam ${bam_loc}GSM1603209.filt.nodup.bam ${bam_loc}GSM1603215.filt.nodup.bam ${bam_loc}GSM1603227.filt.nodup.bam ${bam_loc}GSM1603211.filt.nodup.bam ${bam_loc}GSM1603217.filt.nodup.bam ${bam_loc}GSM1603223.filt.nodup.bam ${bam_loc}GSM1603221.filt.nodup.bam \
  --extendReads 100 \
  --labels SRR2043614_H3K27ac SRR1603646_H3K27ac SRR1603650_H3K27ac SRR1603654_H3K27ac SRR1522115_med1 SRR061743_H3K4m3 SRR061747_H3K79me2 SRR647929_H3K27me3 SRR577482_H3K4me3 SRR577483_H3K4me3 SRR074201_gammaH2AX SRR074195_H2AX SRR074196_H2AX SRR3722566_BRD4 SRR3722577_Pol2_DMSO_4h SRR6010201_YY1 SRR2157609_DMSO_PolII GSM1603225_GFP_H3K4me1 GSM1603213_GFP_H3K4me3 GSM1603209_GFP_H3K36me3 GSM1603215_GFP_H3K79me3 GSM1603227_GFP_H3K9me3 GSM1603211_GFP_H3K27Ac GSM1603217_GFP_PolII GSM1603223_GFP_S5P GSM1603221_GFP_S2P \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/experiments_mbs.npz

plotPCA \
  -in /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/experiments_mbs.npz \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/PCA_readCounts_experiments.pdf \
  -T "PCA of read counts - Experiments"

### H3K27
bam_loc='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/'
multiBamSummary bins \
  --bamfiles ${bam_loc}SRR2043614.filt.nodup.bam ${bam_loc}SRR1603646.filt.nodup.bam ${bam_loc}SRR1603650.filt.nodup.bam ${bam_loc}SRR1603654.filt.nodup.bam ${bam_loc}GSM1603211.filt.nodup.bam ${bam_loc}SRR647929.filt.nodup.bam \
  --extendReads 100 \
  --labels SRR2043614_H3K27ac SRR1603646_H3K27ac SRR1603650_H3K27ac SRR1603654_H3K27ac GSM1603211_GFP_H3K27Ac SRR647929_H3K27me3  \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/H3K27_mbs.npz

plotPCA \
  -in /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/H3K27_mbs.npz \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/PCA_readCounts_H3K27.pdf \
  -T "PCA of read counts - H3K27"

### H3K4
bam_loc='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/'
multiBamSummary bins \
  --bamfiles ${bam_loc}SRR061743.filt.nodup.bam ${bam_loc}SRR577482.filt.nodup.bam ${bam_loc}SRR577483.filt.nodup.bam ${bam_loc}GSM1603225.filt.nodup.bam ${bam_loc}GSM1603213.filt.nodup.bam \
  --extendReads 100 \
  --labels SRR061743_H3K4m3 SRR577482_H3K4me3 SRR577483_H3K4me3 GSM1603225_GFP_H3K4me1 GSM1603213_GFP_H3K4me3 \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/H3K4_mbs.npz

plotPCA \
  -in /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/H3K4_mbs.npz \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/PCA_readCounts_H3K4.pdf \
  -T "PCA of read counts - H3K4"

### PolII
bam_loc='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/'
multiBamSummary bins \
  --bamfiles ${bam_loc}SRR1522115.filt.nodup.bam ${bam_loc}SRR3722566.filt.nodup.bam ${bam_loc}SRR3722577.filt.nodup.bam ${bam_loc}SRR2157609.filt.nodup.bam ${bam_loc}GSM1603217.filt.nodup.bam ${bam_loc}GSM1603223.filt.nodup.bam ${bam_loc}GSM1603221.filt.nodup.bam \
  --extendReads 100 \
  --labels SRR1522115_med1 SRR3722566_BRD4 SRR3722577_Pol2_DMSO_4h SRR2157609_DMSO_PolII GSM1603217_GFP_PolII GSM1603223_GFP_S5P GSM1603221_GFP_S2P \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/PolII_mbs.npz

plotPCA \
  -in /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/PolII_mbs.npz \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/PCA_readCounts_PolII.pdf \
  -T "PCA of read counts - PolII"






