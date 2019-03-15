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

ddsFullCountTable$condition <- relevel(ddsFullCountTable$condition, "damctl")
dds <- DESeq(ddsFullCountTable)
res <- results(dds)
resSig <- res[which(res$padj < 0.05),]

write.table(res, "damIDseq.txt")
write.table(resSig, "damIDseq_sig0.05.txt")
