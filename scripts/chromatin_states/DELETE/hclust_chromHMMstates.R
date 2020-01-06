### This script runs hierarchial clustering on the chromHMM output
library(pheatmap)
setwd("/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/dendogram_for_reorder")
cHMMe <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder/emissions_15.txt", header = T, sep="\t", row.names=1)

### Original
cHMMe <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/emissions_15_labeled.txt", header = T, sep="\t", row.names=1)
d <- dist(cHMMe, method="euclidean")
hc1 <- hclust(d, method="complete")
pdf("chromHMM_dendogram_original.pdf")
plot(hc1)
dev.off()
pdf("chromHMM_heatmap_original.pdf")
pheatmap(cHMMe)
dev.off()

### Reorder 1
cHMMe <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder/emissions_15.txt", header = T, sep="\t", row.names=1)
d <- dist(cHMMe, method="euclidean")
hc1 <- hclust(d, method="complete")
pdf("chromHMM_dendogram_reorder1.pdf")
plot(hc1)
dev.off()
pdf("chromHMM_heatmap_reorder1.pdf")
pheatmap(cHMMe)
dev.off()

### Reorder 2
cHMMe <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/emissions_15_reorder2.txt", header = T, sep="\t", row.names=1)
d <- dist(cHMMe, method="euclidean")
hc1 <- hclust(d, method="complete")
pdf("chromHMM_dendogram_reorder2.pdf")
plot(hc1)
dev.off()
pdf("chromHMM_heatmap_reorder2.pdf")
pheatmap(cHMMe)
dev.off()

### dendogram order; reorder3
cHMMe <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/emissions_15_reorder3_cluster.txt", header = T, sep="\t", row.names=1)
d <- dist(cHMMe, method="euclidean")
hc1 <- hclust(d, method="complete")
pdf("chromHMM_dendogram_cluster_order.pdf")
plot(hc1)
dev.off()
pdf("chromHMM_heatmap_cluster_order.pdf")
pheatmap(cHMMe)
dev.off()
