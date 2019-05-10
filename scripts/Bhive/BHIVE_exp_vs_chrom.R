### This R script makes a violin scatter plot by x=chromosome, y=Hiv expression

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots')

### Load file
hivexp <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt", header=T, sep='\t')

### Remove line with unknown chromosome
hivexp <- hivexp[!hivexp$chr == "chrUn_GL000195v1", ]

### Make plot
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(scales)

hivexp$chr <- factor(hivexp$chr,levels = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                           "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                                           "chr21","chr22","chrX","chrY"))
pdf("BHIVE_exp_by_chrom.pdf")
ggplot(hivexp, aes(x=chr, y=expr)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
  ggtitle("HIV Expression by Chromosome") + 
  ylab("HIV Expression in log10") +
  #stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.5, color="red") +
  scale_x_discrete(labels = wrap_format(10)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()