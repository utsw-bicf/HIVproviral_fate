### This script runs the differential analysis of deseq2 on the damIDseq following
### https://github.com/sarahhcarl/flychip/wiki/Basic-DamID-analysis-pipeline; Differential enrichment analysis
library("DESeq2")

setwd("/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/DE")
### Read file
df <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/DE/combine_GATCcount_table.tsv", header=T, sep="\t")
ct <- df[,2:length(df)]

DamID_design <- data.frame(row.names = c("LaminB11", "LaminB12", "Dam1", "Dam2"),
                           condition=c("LaminB1", "LaminB1", "damctl", "damctl"),
                           libType = c("single-end", "single-end", "single-end", "single-end"))


ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = ct,
  colData = DamID_design,
  design = ~ condition)