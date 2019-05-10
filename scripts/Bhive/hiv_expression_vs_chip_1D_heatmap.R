### This script makes a heatmap of HIV expression of 5 groups
### This will be placed next to the chip_1D graphics

library(pheatmap)
library(RColorBrewer)

setwd("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip")

### Load the files
g1 <- read.table("number1.bed", header=F, sep="\t")
g2 <- read.table("number2.bed", header=F, sep="\t")
g3 <- read.table("number3.bed", header=F, sep="\t")
g4 <- read.table("number4.bed", header=F, sep="\t")
g5 <- read.table("number5.bed", header=F, sep="\t")

### Merge and reduce files
df <- rbind(g4, g5, g1, g2, g3)
dfr <- df$V5

### Make matrix
mat <- as.matrix(dfr)

### Make heatmap
pheatmap(mat, filename = "BHIVE_expression_group_heatmap.pdf", 
         cluster_cols = F, 
         cluster_rows = F, 
         gaps_row = c(581, 1223, 1350, 1421),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cellwidth = 20)
dev.off()

