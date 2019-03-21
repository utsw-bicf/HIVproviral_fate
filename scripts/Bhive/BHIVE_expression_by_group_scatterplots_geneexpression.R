### this R script makes the scatter plots for bhive results
library(ggplot2)
library(tidyr)
library(RColorBrewer)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots')

########################################
### Load bed files by group
g1 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed", header =F, sep="\t")
g2 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed", header =F, sep="\t")
g3 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed", header =F, sep="\t")
g4 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed", header =F, sep="\t")
g5 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed", header =F, sep="\t")

### add groups and concatenate tables
g1$group <- "Intergenic_Same"
g2$group <- "Intergenic_Convergent"
g3$group <- "Intergenic_Divergent"
g4$group <- "Intragenic_Same"
g5$group <- "Intragenic_Convergent"

### Make correct ensembl number
g1$ENSEMBL <- substr(g1$V15, 0, 15)
g2$ENSEMBL <- substr(g2$V15, 0, 15)
g3$ENSEMBL <- substr(g3$V15, 0, 15)
g4$ENSEMBL <- substr(g4$V15, 0, 15)
g5$ENSEMBL <- substr(g5$V15, 0, 15)

### Load RNAseq data
logCPM <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq_reg/workflow/output_PE/countTable.logCPM.txt", header=T, sep="\t")
logCPM$mean <- rowMeans(logCPM[c('Emily_jurkat_A', 'Emily_jurkat_B', 'Emily_jurkat_C')], na.rm=TRUE)

fpkm <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq_reg/workflow/output_PE/countTable.fpkm.txt", header=T, sep="\t")
fpkm$log10fpkm <- log(rowMeans(fpkm[c('Emily_jurkat_A', 'Emily_jurkat_B', 'Emily_jurkat_C')], na.rm=TRUE)+1)

### Reduce rna seq and add to groups table
rfpkm <- fpkm[,c("ENSEMBL", "log10fpkm")]
g1fpkm <- merge(g1, rfpkm, by="ENSEMBL")
g2fpkm <- merge(g2, rfpkm, by="ENSEMBL")
g3fpkm <- merge(g3, rfpkm, by="ENSEMBL")
g4fpkm <- merge(g4, rfpkm, by="ENSEMBL")
g5fpkm <- merge(g5, rfpkm, by="ENSEMBL")

### Pearson's correlation
cor_g1fpkm <- cor.test(g1fpkm$log10fpkm, g1fpkm$V11, method = "pearson", conf.level = 0.95)
cor_g2fpkm <- cor.test(g2fpkm$log10fpkm, g2fpkm$V11, method = "pearson", conf.level = 0.95)
cor_g3fpkm <- cor.test(g3fpkm$log10fpkm, g3fpkm$V11, method = "pearson", conf.level = 0.95)
cor_g4fpkm <- cor.test(g4fpkm$log10fpkm, g4fpkm$V11, method = "pearson", conf.level = 0.95)
cor_g5fpkm <- cor.test(g5fpkm$log10fpkm, g5fpkm$V11, method = "pearson", conf.level = 0.95)

### Merge Tables
dffpkm <- rbind(g1fpkm, g2fpkm, g3fpkm, g4fpkm, g5fpkm)



### Make scatter plot
### add pearson correlation manually
### "#CC79A7", "#D55E00", "#0072B2", "#F0E442", "#009E73"
pdf("BHIVE_expression_scatterplot_expression_all.pdf")
ggplot(dffpkm, aes(x=V11, y=log10fpkm, size=V8, color=group)) + 
  scale_color_manual(values=c("#CC79A7", "#D55E00", "#0072B2", "#F0E442", "#009E73")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Host Gene Expression") + 
  xlab("Log10 HIV Expression") + 
  ylab("Host Gene Expression (Log10(fpkm))") + 
  scale_x_continuous(limits=c(-1.5, 1.5)) +
  scale_y_continuous(limits=c(0, 6.5)) +
  annotate("text", x = -1.5, y = 6.40, label = "Intergenic Same r = 0.18 (p-value = 0.048)", size=3, color = "#0072B2", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = -1.5, y = 6.25, label = "Intergenic Convergent r = 0.12 (p-value = 0.318)", size=3, color = "#CC79A7", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = -1.5, y = 6.10, label = "Intergenic Divergent r = 0.11 (p-value = 0.326)", size=3, color = "#D55E00", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = -1.5, y = 5.95, label = "Intragenic Same r = 0.02 (p-value = 0.601)", size=3, color = "#009E73", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = -1.5, y = 5.8, label = "Intragenic Convergent r = 0.04 (p-value = 0.113)", size=3, color = "#F0E442", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

##################################################################################
### Make individual scatterplots
pdf("BHIVE_expression_scatterplot_expression_group1.pdf")
ggplot(g1fpkm, aes(x=V11, y=log10fpkm, size=V8, color=group)) + 
  scale_color_manual(values=c("#0072B2")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Host Gene Expression") + 
  xlab("Log10 HIV Expression") + 
  ylab("Host Gene Expression (Log10(fpkm))") + 
  scale_x_continuous(limits=c(-1.5, 1.5)) +
  scale_y_continuous(limits=c(0, 6.5)) +
  annotate("text", x = -1.5, y = 6.40, label = "Intergenic Same r = 0.18 (p-value = 0.048)", size=3, color = "#0072B2", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_expression_group2.pdf")
ggplot(g2fpkm, aes(x=V11, y=log10fpkm, size=V8, color=group)) + 
  scale_color_manual(values=c("#CC79A7")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Host Gene Expression") + 
  xlab("Log10 HIV Expression") + 
  ylab("Host Gene Expression (Log10(fpkm))") + 
  scale_x_continuous(limits=c(-1.5, 1.5)) +
  scale_y_continuous(limits=c(0, 6.5)) +
  annotate("text", x = -1.5, y = 6.40, label = "Intergenic Convergent r = 0.12 (p-value = 0.318)", size=3, color = "#CC79A7", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_expression_group3.pdf")
ggplot(g3fpkm, aes(x=V11, y=log10fpkm, size=V8, color=group)) + 
  scale_color_manual(values=c("#D55E00")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Host Gene Expression") + 
  xlab("Log10 HIV Expression") + 
  ylab("Host Gene Expression (Log10(fpkm))") + 
  scale_x_continuous(limits=c(-1.5, 1.5)) +
  scale_y_continuous(limits=c(0, 6.5)) +
  annotate("text", x = -1.5, y = 6.40, label = "Intergenic Divergent r = 0.11 (p-value = 0.326)", size=3, color = "#D55E00", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_expression_group4.pdf")
ggplot(g4fpkm, aes(x=V11, y=log10fpkm, size=V8, color=group)) + 
  scale_color_manual(values=c("#009E73")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Host Gene Expression") + 
  xlab("Log10 HIV Expression") + 
  ylab("Host Gene Expression (Log10(fpkm))") + 
  scale_x_continuous(limits=c(-1.5, 1.5)) +
  scale_y_continuous(limits=c(0, 6.5)) +
  annotate("text", x = -1.5, y = 6.40, label = "Intragenic Same r = 0.02 (p-value = 0.601)", size=3, color = "#009E73", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_expression_group5.pdf")
ggplot(g5fpkm, aes(x=V11, y=log10fpkm, size=V8, color=group)) + 
  scale_color_manual(values=c("#F0E442")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Host Gene Expression") + 
  xlab("Log10 HIV Expression") + 
  ylab("Host Gene Expression (Log10(fpkm))") + 
  scale_x_continuous(limits=c(-1.5, 1.5)) +
  scale_y_continuous(limits=c(0, 6.5)) +
  annotate("text", x = -1.5, y = 6.40, label = "Intragenic Convergent r = 0.04 (p-value = 0.113)", size=3, color = "#F0E442", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()



