### This script runs hierarchial clustering on the chromHMM output
library(pheatmap)
setwd("/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15")

### Original
cHMMe <- read.table(file="emissions_15.txt", header = T, sep="\t", row.names=1)
d <- dist(cHMMe, method="euclidean")
hc1 <- hclust(d, method="complete")
pdf("chromHMM_dendogram_original.pdf")
plot(hc1)
dev.off()
pdf("chromHMM_heatmap_original.pdf")
pheatmap(cHMMe)
dev.off()

