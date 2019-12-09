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

### add log10(distance)
g1$logDis <- log10(g1$V23)
g2$logDis <- log10(g2$V23)
g3$logDis <- log10(g3$V23)
g4$logDis <- log10(g4$V23)
g5$logDis <- log10(g5$V23)

### add groups
g1$group <- "Intergenic_Same"
g2$group <- "Intergenic_Convergent"
g3$group <- "Intergenic_Divergent"
g4$group <- "Intragenic_Same"
g5$group <- "Intragenic_Convergent"

### Pearson's correlation
cor_g1 <- cor.test(g1$logDis, g1$V11, method = "pearson", conf.level = 0.95)
cor_g2 <- cor.test(g2$logDis, g2$V11, method = "pearson", conf.level = 0.95)
cor_g3 <- cor.test(g3$logDis, g3$V11, method = "pearson", conf.level = 0.95)
cor_g4 <- cor.test(g4$logDis, g4$V11, method = "pearson", conf.level = 0.95)
cor_g5 <- cor.test(g5$logDis, g5$V11, method = "pearson", conf.level = 0.95)

### Merge Tables
df <- rbind(g1, g2, g3, g4, g5)

### Make scatter plot
### add pearson correlation manually
### "#ff0000", "#D55E00", "#0072B2", "#808080", "#009E73"
pdf("BHIVE_expression_scatterplot_closestTSS_all.pdf")
ggplot(df, aes(x=logDis, y=V11, size=V8, color=group)) + 
  scale_color_manual(values=c("#ff0000", "#D55E00", "#0072B2", "#808080", "#009E73")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("Correlation of HIV Expression to Nearest TSS") + 
  xlab("Log10 Distance to TSS") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(2, 6.5)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 2.0, y = 3, label = "Intergenic Same r = -0.13 (p-value = 0.14)", size=3, color = "#0072B2", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = 2.0, y = 2.85, label = "Intergenic Convergent r = -0.19 (p-value = 0.11)", size=3, color = "#ff0000", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = 2.0, y = 2.70, label = "Intergenic Divergent r = -0.16 (p-value = 0.15)", size=3, color = "#D55E00", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = 2.0, y = 2.55, label = "Intragenic Same r = -0.13 (p-value = 0.0012)", size=3, color = "#009E73", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = 2.0, y = 2.40, label = "Intragenic Convergent r = -0.23 (p-value = 0.0000000064)", size=3, color = "#808080", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

##################################################################################
### Make individual scatterplots
pdf("BHIVE_expression_scatterplot_closestTSS_group1.pdf")
ggplot(g1, aes(x=logDis, y=V11, size=V8, color=group)) + 
  scale_color_manual(values=c("#0072B2")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("Correlation of HIV Expression to Nearest TSS") + 
  xlab("Log10 Distance to TSS") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(2, 6.5)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 2.0, y = 3, label = "Intergenic Same r = -0.13 (p-value = 0.14)", size=3, color = "#0072B2", hjust=0) + ###This line needs to be manually changed
    theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_closestTSS_group2.pdf")
ggplot(g2, aes(x=logDis, y=V11, size=V8, color=group)) + 
  scale_color_manual(values=c("#ff0000")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("Correlation of HIV Expression to Nearest TSS") + 
  xlab("Log10 Distance to TSS") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(2, 6.5)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 2.0, y = 3, label = "Intergenic Convergent r = -0.19 (p-value = 0.11)", size=3, color = "#ff0000", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_closestTSS_group3.pdf")
ggplot(g3, aes(x=logDis, y=V11, size=V8, color=group)) + 
  scale_color_manual(values=c("#D55E00")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("Correlation of HIV Expression to Nearest TSS") + 
  xlab("Log10 Distance to TSS") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(2, 6.5)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 2.0, y = 3, label = "Intergenic Divergent r = -0.16 (p-value = 0.15)", size=3, color = "#D55E00", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_closestTSS_group4.pdf")
ggplot(g4, aes(x=logDis, y=V11, size=V8, color=group)) + 
  scale_color_manual(values=c("#009E73")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("Correlation of HIV Expression to Nearest TSS") + 
  xlab("Log10 Distance to TSS") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(2, 6.5)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 2.0, y = 3, label = "Intragenic Same r = -0.13 (p-value = 0.0012)", size=3, color = "#009E73", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_closestTSS_group5.pdf")
ggplot(g5, aes(x=logDis, y=V11, size=V8, color=group)) + 
  scale_color_manual(values=c("#808080")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("Correlation of HIV Expression to Nearest TSS") + 
  xlab("Log10 Distance to TSS") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(2, 6.5)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 2.0, y = 3, label = "Intragenic Convergent r = -0.23 (p-value = 0.0000000064)", size=3, color = "#808080", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()
