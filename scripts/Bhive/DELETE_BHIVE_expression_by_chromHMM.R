### This R script makes a violin plot
### X = chromHMM state
### Y = HIV expression

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_chromstates')


### Load BHIVE expression
hivexp <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt", header=T, sep='\t')

### Load chromHMM; 15 states, histone only
chromHMM <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/jurkat_15_dense_chr.bed", header=F, sep='\t', skip=1)

### Add state to BHIVE expression
library(IRanges)
ihivexp <- with(hivexp, IRanges(locus, width=1, names=chr))
ichromHMMs <- with(chromHMM, IRanges(V2, V3, names=V1))
olaps <- findOverlaps(ihivexp, ichromHMMs)
df <- cbind(hivexp[queryHits(olaps),], chromHMM[subjectHits(olaps),])
df2 <- test[which(as.character(test$chr) == as.character(test$V1)),]

### Plot
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(scales)

df2$V4 <- factor(df2$V4,levels = c("1", "2","3","4","5","6","7","8","9","10", "11","12","13","14","15"))

pdf("BHIVE_exp_by_chromstates_histone.pdf")
ggplot(df2, aes(x=V4, y=expr)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0)) +
  ggtitle("chromHMM state of HIV expression") + 
  ylab("HIV Expression") +
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.5, color="red") +
  scale_x_discrete(labels = wrap_format(10)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

######################################################################
######################################################################
######################################################################

### Do the same as above but with all markers
### Load BHIVE expression
hivexp <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt", header=T, sep='\t')

### Load chromHMM; 15 states, histone only
chromHMM <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/states_all_data/allbam_learn15/jurkat_15_dense_chr.bed", header=F, sep='\t', skip=1)

### Add state to BHIVE expression
library(IRanges)
ihivexp <- with(hivexp, IRanges(locus, width=1, names=chr))
ichromHMMs <- with(chromHMM, IRanges(V2, V3, names=V1))
olaps <- findOverlaps(ihivexp, ichromHMMs)
df <- cbind(hivexp[queryHits(olaps),], chromHMM[subjectHits(olaps),])
df2 <- test[which(as.character(test$chr) == as.character(test$V1)),]

### Plot
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(scales)

df2$V4 <- factor(df2$V4,levels = c("1", "2","3","4","5","6","7","8","9","10", "11","12","13","14","15"))

pdf("BHIVE_exp_by_chromstates_all.pdf")
ggplot(df2, aes(x=V4, y=expr)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0)) +
  ggtitle("chromHMM state of HIV expression") + 
  ylab("HIV Expression") +
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.5, color="red") +
  scale_x_discrete(labels = wrap_format(10)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()