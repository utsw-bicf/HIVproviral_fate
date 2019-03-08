### this R script makes the scatter plots for bhive results
library(ggplot2)
#require(gridExtra)
library(tidyr)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate')

################################################### Ivan's Scatter Plots ###################################################
### Load the same and opposite bed files
same <- read.table("hiv_expression_proximalTSS_same_plus_dis.bed", header=F, sep="\t")
opp <- read.table("hiv_expression_proximalTSS_opposite_plus_dis.bed", header=F, sep="\t")

### add log10(distance)
same$logDis <- log10(same$V23)
opp$logDis <- log10(opp$V23)

### Make scatter plot
pdf("BHIVE_expression_scatterplot_same.pdf")
ggplot(same, aes(x=logDis, y=V11, size=V8)) + 
  geom_point(alpha=1/4) + 
  ggtitle("Same Orientation") + 
  xlab("Log10 Distance to TSS") + 
  ylab("HIV Expression") + 
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()


pdf("BHIVE_expression_scatterplot_opposite.pdf")
ggplot(opp, aes(x=logDis, y=V11, size=V8)) + 
  geom_point(alpha=1/4, color="red") + 
  ggtitle("Opposite Orientation") + 
  xlab("Log10 Distance to TSS") + 
  ylab("HIV Expression") + 
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()



################################################### Genic Scatter Plots ###################################################
### Load the same and opposite bed files
same1 <- read.table("hiv_expression_genic_same_plus_dis.bed", header=F, sep="\t")
opp1 <- read.table("hiv_expression_genic_opposite_plus_dis.bed", header=F, sep="\t")

### Split the ensembl ID name
same2 <- same1 %>% separate(V15, c("ENSEMBL", "B"))
opp2 <- opp1 %>% separate(V15, c("ENSEMBL", "B"))

### add log10(distance)
same2$logDis <- log10(same2$V22)
opp2$logDis <- log10(opp2$V22)

### Load RNAseq data
logCPM <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq_reg/workflow/output_PE/countTable.logCPM.txt", header=T, sep="\t")
logCPM$avg3 <- rowMeans(logCPM[,c("Emily_jurkat_A", "Emily_jurkat_B", "Emily_jurkat_C")])

### Merge data frames
same3 <- merge(same2, logCPM, by="ENSEMBL", all.x=T)
opp3 <- merge(opp2, logCPM, by="ENSEMBL", all.x=T)


### Make scatter plot
pdf("BHIVE_expression_scatterplot_same_genic.pdf")
ggplot(same3, aes(x=V11, y=avg3, size=logDis, na.rm = TRUE)) + 
  geom_point(alpha=1/4) + 
  ggtitle("Same Orientation") + 
  xlab("HIV Expression") + 
  ylab("RNA-seq logCPM") + 
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()


pdf("BHIVE_expression_scatterplot_opposite_genic.pdf")
ggplot(opp3, aes(x=V11, y=avg3, size=logDis, na.rm = TRUE)) + 
  geom_point(alpha=1/4) + 
  ggtitle("Opposite Orientation") + 
  xlab("HIV Expression") + 
  ylab("RNA-seq logCPM") + 
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

################################################### Exonic Scatter Plots ###################################################
### Load the same and opposite bed files
esame1 <- read.table("hiv_expression_exon_same_plus_dis.bed", header=F, sep="\t")
eopp1 <- read.table("hiv_expression_exon_opposite_plus_dis.bed", header=F, sep="\t")

### Split the ensembl ID name
esame2 <- esame1 %>% separate(V15, c("ENSEMBL", "B"))
eopp2 <- eopp1 %>% separate(V15, c("ENSEMBL", "B"))

### add log10(distance)
esame2$logDis <- log10(esame2$V22)
eopp2$logDis <- log10(eopp2$V22)

### Load RNAseq data
logCPM <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq_reg/workflow/output_PE/countTable.logCPM.txt", header=T, sep="\t")
logCPM$avg3 <- rowMeans(logCPM[,c("Emily_jurkat_A", "Emily_jurkat_B", "Emily_jurkat_C")])

### Merge data frames
esame3 <- merge(esame2, logCPM, by="ENSEMBL", all.x=T)
eopp3 <- merge(eopp2, logCPM, by="ENSEMBL", all.x=T)


### Make scatter plot
pdf("BHIVE_expression_scatterplot_same_exon.pdf")
ggplot(esame3, aes(x=V11, y=avg3, size=logDis, na.rm = TRUE)) + 
  geom_point(alpha=1/4) + 
  ggtitle("Same Orientation") + 
  xlab("HIV Expression") + 
  ylab("RNA-seq logCPM") + 
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()


pdf("BHIVE_expression_scatterplot_opposite_exon.pdf")
ggplot(eopp3, aes(x=V11, y=avg3, size=logDis, na.rm = TRUE)) + 
  geom_point(alpha=1/4) + 
  ggtitle("Opposite Orientation") + 
  xlab("HIV Expression") + 
  ylab("RNA-seq logCPM") + 
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()
