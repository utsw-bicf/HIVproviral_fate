library(ggplot2)
library(tidyr)
library(RColorBrewer)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestEnhancers_scatterplots')

########################################
### Load bed files by group
g1 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same_closestEnhancer.bed", header =F, sep="\t")
g2 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite_closestEnhancer.bed", header =F, sep="\t")
g3 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite_closestEnhancer.bed", header =F, sep="\t")
g4 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene_closestEnhancer.bed", header =F, sep="\t")
g5 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene_closestEnhancer.bed", header =F, sep="\t")
g6 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6thru8_2gene_closestEnhancer.bed", header =F, sep="\t")

### add log10(distance)
g1$logDis <- log10(abs(g1$V27))
g2$logDis <- log10(abs(g2$V27))
g3$logDis <- log10(abs(g3$V27))
g4$logDis <- log10(abs(g4$V27))
g5$logDis <- log10(abs(g5$V27))
g6$logDis <- log10(abs(g6$V20))

### add groups
g1$group <- "Intergenic_Same"
g2$group <- "Intergenic_Convergent"
g3$group <- "Intergenic_Divergent"
g4$group <- "Intragenic_Same"
g5$group <- "Intragenic_Convergent"
g6$group <- "Overlaps Two Genes"

### Pearson's correlation
cor_g1 <- cor.test(g1$logDis, g1$V11, method = "pearson", conf.level = 0.95)
cor_g2 <- cor.test(g2$logDis, g2$V11, method = "pearson", conf.level = 0.95)
cor_g3 <- cor.test(g3$logDis, g3$V11, method = "pearson", conf.level = 0.95)
cor_g4 <- cor.test(g4$logDis, g4$V11, method = "pearson", conf.level = 0.95)
cor_g5 <- cor.test(g5$logDis, g5$V11, method = "pearson", conf.level = 0.95)
cor_g6 <- cor.test(g6$logDis, g6$V15, method = "pearson", conf.level = 0.95)

### Merge Tables
g1r <- subset(g1,select=c(V11,logDis,V8,group))
colnames(g1r) <- c("exp", "logDis", "mapq", "group")
g2r <- subset(g2,select=c(V11,logDis,V8,group))
colnames(g2r) <- c("exp", "logDis", "mapq", "group")
g3r <- subset(g3,select=c(V11,logDis,V8,group))
colnames(g3r) <- c("exp", "logDis", "mapq", "group")
g4r <- subset(g4,select=c(V11,logDis,V8,group))
colnames(g4r) <- c("exp", "logDis", "mapq", "group")
g5r <- subset(g5,select=c(V11,logDis,V8,group))
colnames(g5r) <- c("exp", "logDis", "mapq", "group")
g6r <- subset(g6,select=c(V15,logDis,V12,group))
colnames(g6r) <- c("exp", "logDis", "mapq", "group")
df <- rbind(g1r, g2r, g3r, g4r, g5r,g6r)

### Make scatter plot
### add pearson correlation manually
### "#CC79A7", "#D55E00", "#0072B2", "#F0E442", "#009E73"
pdf("BHIVE_expression_scatterplot_closestEnhancer_all.pdf")
ggplot(df, aes(x=logDis, y=exp, size=mapq, color=group)) + 
  scale_color_manual(values=c("#CC79A7", "#D55E00", "#0072B2", "#F0E442", "#009E73","#9933CC")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("Correlation of HIV Expression to Nearest Enhancer") + 
  xlab("Log10 Distance to Enhancer") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(2.5, 7)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 2.5, y = 3, label = "Intergenic Same r = -0.13 (p-value = 0.14)", size=3, color = "#0072B2", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = 2.5, y = 2.85, label = "Intergenic Convergent r = -0.05 (p-value = 0.69)", size=3, color = "#CC79A7", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = 2.5, y = 2.70, label = "Intergenic Divergent r = -0.12 (p-value = 0.32)", size=3, color = "#D55E00", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = 2.5, y = 2.55, label = "Intragenic Same r = -0.01 (p-value = 0.80)", size=3, color = "#009E73", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = 2.5, y = 2.40, label = "Intragenic Convergent r = -0.03 (p-value = 0.38)", size=3, color = "#F0E442", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = 2.5, y = 2.25, label = "Overlaps 2 Genes r = 0.14 (p-value = 0.27)", size=3, color = "#9933CC", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

##################################################################################
### Make individual scatterplots
pdf("BHIVE_expression_scatterplot_closestEnhancer_group1.pdf")
ggplot(g1r, aes(x=logDis, y=exp, size=mapq, color=group)) + 
  scale_color_manual(values=c("#0072B2")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("Correlation of HIV Expression to Nearest Enhancer") + 
  xlab("Log10 Distance to Enhancer") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(2.5, 7)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 2.5, y = 3, label = "Intergenic Same r = -0.13 (p-value = 0.14)", size=3, color = "#0072B2", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_closestEnhancer_group2.pdf")
ggplot(g2r, aes(x=logDis, y=exp, size=mapq, color=group)) + 
  scale_color_manual(values=c("#CC79A7")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("Correlation of HIV Expression to Nearest Enhancer") + 
  xlab("Log10 Distance to Enhancer") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(2.5, 7)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 2.5, y = 3, label = "Intergenic Convergent r = -0.05 (p-value = 0.69)", size=3, color = "#CC79A7", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_closestEnhancer_group3.pdf")
ggplot(g3r, aes(x=logDis, y=exp, size=mapq, color=group)) + 
  scale_color_manual(values=c("#D55E00")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("Correlation of HIV Expression to Nearest Enhancer") + 
  xlab("Log10 Distance to Enhancer") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(2.5, 7)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 2.5, y = 3, label = "Intergenic Divergent r = -0.12 (p-value = 0.32)", size=3, color = "#D55E00", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_closestEnhancer_group4.pdf")
ggplot(g4r, aes(x=logDis, y=exp, size=mapq, color=group)) + 
  scale_color_manual(values=c("#009E73")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("Correlation of HIV Expression to Nearest Enhancer") + 
  xlab("Log10 Distance to Enhancer") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(2.5, 7)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 2.5, y = 3, label = "Intragenic Same r = -0.01 (p-value = 0.80)", size=3, color = "#009E73", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_closestEnhancer_group5.pdf")
ggplot(g5r, aes(x=logDis, y=exp, size=mapq, color=group)) + 
  scale_color_manual(values=c("#F0E442")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("Correlation of HIV Expression to Nearest Enhancer") + 
  xlab("Log10 Distance to Enhancer") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(2.5, 7)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 2.5, y = 3, label = "Intragenic Convergent r = -0.03 (p-value = 0.38)", size=3, color = "#F0E442", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_closestEnhancer_group6.pdf")
ggplot(g6r, aes(x=logDis, y=exp, size=mapq, color=group)) + 
  scale_color_manual(values=c("#9933CC")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("Correlation of HIV Expression to Nearest Enhancer") + 
  xlab("Log10 Distance to Enhancer") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(2.5, 7)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 2.5, y = 3, label = "Overlaps 2 Genes r = 0.14 (p-value = 0.27)", size=3, color = "#9933CC", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()
