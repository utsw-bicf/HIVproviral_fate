library(ggplot2)
library(tidyr)
library(RColorBrewer)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestEnhancers_scatterplots')

########################################
### Load bed files by group
g1 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same_closestSuperEnhancer.bed", header =F, sep="\t")
g2 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite_closestSuperEnhancer.bed", header =F, sep="\t")
g3 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite_closestSuperEnhancer.bed", header =F, sep="\t")
g4 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene_closestSuperEnhancer.bed", header =F, sep="\t")
g5 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene_closestSuperEnhancer.bed", header =F, sep="\t")
g6 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6thru8_2gene_closestSuperEnhancer.bed", header =F, sep="\t")

### add log10(distance)
g1$logDis <- log10(abs(g1$V27))
is.na(g1)<-sapply(g1, is.infinite)
g1[is.na(g1)]<-0
g1 <- g1[!grepl("chrY",g1$V1),]

g2$logDis <- log10(abs(g2$V27))
is.na(g2)<-sapply(g2, is.infinite)
g2[is.na(g2)]<-0
g2 <- g2[!grepl("chrY",g2$V1),]

g3$logDis <- log10(abs(g3$V27))
is.na(g3)<-sapply(g3, is.infinite)
g3[is.na(g3)]<-0
g3 <- g3[!grepl("chrY",g3$V1),]

g4$logDis <- log10(abs(g4$V27))
is.na(g4)<-sapply(g4, is.infinite)
g4[is.na(g4)]<-0
g4 <- g4[!grepl("chrY",g4$V1),]

g5$logDis <- log10(abs(g5$V27))
is.na(g5)<-sapply(g5, is.infinite)
g5[is.na(g5)]<-0
g5 <- g5[!grepl("chrY",g5$V1),]

g6$logDis <- log10(abs(g6$V20))
is.na(g6)<-sapply(g6, is.infinite)
g6[is.na(g6)]<-0
g6 <- g6[!grepl("chrY",g6$V1),]

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
### "#ff0000", "#D55E00", "#0072B2", "#808080", "#009E73"
pdf("BHIVE_expression_scatterplot_closestSuperEnhancer_all.pdf")
ggplot(df, aes(x=logDis, y=exp, size=mapq, color=group)) + 
  scale_color_manual(values=c("#ff0000", "#D55E00", "#0072B2", "#808080", "#009E73","#9933CC")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Nearest SuperEnhancer") + 
  xlab("Log10 Distance to SuperEnhancer") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(0, 7.5)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 0, y = 3, label = "Intergenic Same r = -0.02 (p-value = 0.87)", size=3, color = "#0072B2", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = 0, y = 2.85, label = "Intergenic Convergent r = -0.08 (p-value = 0.51)", size=3, color = "#ff0000", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = 0, y = 2.70, label = "Intergenic Divergent r = -0.11 (p-value = 0.36)", size=3, color = "#D55E00", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = 0, y = 2.55, label = "Intragenic Same r = -0.09 (p-value = 0.02)", size=3, color = "#009E73", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = 0, y = 2.40, label = "Intragenic Convergent r = -0.12 (p-value = 0.002)", size=3, color = "#808080", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = 0, y = 2.25, label = "Overlaps 2 Genes r = -0.05 (p-value = 0.73)", size=3, color = "#9933CC", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

##################################################################################
### Make individual scatterplots
pdf("BHIVE_expression_scatterplot_closestSuperEnhancer_group1.pdf")
ggplot(g1r, aes(x=logDis, y=exp, size=mapq, color=group)) + 
  scale_color_manual(values=c("#0072B2")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Nearest SuperEnhancer") + 
  xlab("Log10 Distance to SuperEnhancer") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(0, 7.5)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 0, y = 3, label = "Intergenic Same r = -0.02 (p-value = 0.87)", size=3, color = "#0072B2", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_closestSuperEnhancer_group2.pdf")
ggplot(g2r, aes(x=logDis, y=exp, size=mapq, color=group)) + 
  scale_color_manual(values=c("#ff0000")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Nearest SuperEnhancer") + 
  xlab("Log10 Distance to SuperEnhancer") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(0, 7.5)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 0, y = 3, label = "Intergenic Convergent r = -0.08 (p-value = 0.51)", size=3, color = "#ff0000", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_closestSuperEnhancer_group3.pdf")
ggplot(g3r, aes(x=logDis, y=exp, size=mapq, color=group)) + 
  scale_color_manual(values=c("#D55E00")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Nearest SuperEnhancer") + 
  xlab("Log10 Distance to SuperEnhancer") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(0, 7.5)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 0, y = 3, label = "Intergenic Divergent r = -0.11 (p-value = 0.36)", size=3, color = "#D55E00", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_closestSuperEnhancer_group4.pdf")
ggplot(g4r, aes(x=logDis, y=exp, size=mapq, color=group)) + 
  scale_color_manual(values=c("#009E73")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Nearest SuperEnhancer") + 
  xlab("Log10 Distance to SuperEnhancer") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(0, 7.5)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 0, y = 3, label = "Intragenic Same r = -0.09 (p-value = 0.02)", size=3, color = "#009E73", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_closestSuperEnhancer_group5.pdf")
ggplot(g5r, aes(x=logDis, y=exp, size=mapq, color=group)) + 
  scale_color_manual(values=c("#808080")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Nearest SuperEnhancer") + 
  xlab("Log10 Distance to SuperEnhancer") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(0, 7.5)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 0, y = 3, label = "Intragenic Convergent r = -0.12 (p-value = 0.002)", size=3, color = "#808080", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

pdf("BHIVE_expression_scatterplot_closestSuperEnhancer_group6.pdf")
ggplot(g6r, aes(x=logDis, y=exp, size=mapq, color=group)) + 
  scale_color_manual(values=c("#9933CC")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Nearest SuperEnhancer") + 
  xlab("Log10 Distance to SuperEnhancer") + 
  ylab("HIV Expression") + 
  scale_x_continuous(limits=c(0, 7.5)) +
  scale_y_continuous(limits=c(-3.5, 3)) +
  annotate("text", x = 0, y = 3, label = "Overlaps 2 Genes r = -0.05 (p-value = 0.73)", size=3, color = "#9933CC", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()
