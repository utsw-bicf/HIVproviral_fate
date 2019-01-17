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
