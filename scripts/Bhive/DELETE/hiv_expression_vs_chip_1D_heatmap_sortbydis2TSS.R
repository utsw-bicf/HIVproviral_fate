### This creates an expression heatmap sorted by distance to TSS
library(pheatmap)
library(RColorBrewer)

setwd("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip")

### Load the files
g1 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed", header=F, sep="\t")
g2 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed", header=F, sep="\t")
g3 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed", header=F, sep="\t")
g4 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed", header=F, sep="\t")
g5 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed", header=F, sep="\t")

### Sort by distance (V23)
sg1 <- g1[order(g1$V23),]
sg2 <- g2[order(g2$V23),]
sg3 <- g3[order(g3$V23),]
sg4 <- g4[order(g4$V23),]
sg5 <- g5[order(g5$V23),]

### Merge and reduce files
df <- rbind(sg4, sg5, sg1, sg2, sg3)
dfr <- df$V11

### Make matrix
mat <- as.matrix(dfr)

### Make heatmap
pheatmap(mat, filename = "BHIVE_expression_group_heatmap_bydis2TSS.pdf", 
         cluster_cols = F, 
         cluster_rows = F, 
         gaps_row = c(586, 1233, 1361, 1432),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cellwidth = 20)
dev.off()

