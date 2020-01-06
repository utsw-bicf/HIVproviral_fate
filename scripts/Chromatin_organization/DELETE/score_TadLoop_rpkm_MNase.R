#### This script plots the rpkm values of MNase vs HIV expression

library(ggplot2)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_TadLoop')

########## Loops
inloop <- read.table("HIVexp_inLoopScored.txt", header=F, sep="\t")
colnames(inloop) <- c("chrHIV","HIVstart","HIVend", "HIVbrcd", "HIVexp","strand", "chrom",	"start",	"end",	"H3K27ac",	"H3K36me3",	"H3K79me3",	"H3K9me3",	"H3K4me1",	"H3K4me3",	"RNAse",	"MNase",	"DNase",	"TTseq",	"H3K27me3")

### Separate by unique loops
uInLoop <- inloop[!duplicated(inloop[,c("chrHIV","HIVstart","HIVend")]),]
nuInLoop <- inloop[duplicated(inloop[,c("chrHIV","HIVstart","HIVend")]),]

### Plot unique by x = HIV expression y= MNase rpkm
pdf("plots/BHIVE_expression_1x_inLoopscored.pdf")
ggplot(uInLoop, aes(x=HIVexp, y=H3K27me3)) + 
  scale_color_manual(values=c("#0072B2")) +
  geom_point() + 
  #geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Host Gene Expression") + 
  xlab("Log10 HIV Expression") + 
  ylab("MNase-seq rpkm in loop") + 
  scale_x_continuous(limits=c(-3.5, 3)) +
 # scale_y_continuous(limits=c(0, 1)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

### Plot unique by x = HIV expression y= MNase rpkm
pdf("plots/BHIVE_expression_2plus_inLoopscored.pdf")
ggplot(nuInLoop, aes(x=HIVexp, y=H3K27me3)) + 
  scale_color_manual(values=c("#0072B2")) +
  geom_point() + 
  #geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Host Gene Expression") + 
  xlab("Log10 HIV Expression") + 
  ylab("MNase-seq rpkm in loop") + 
  scale_x_continuous(limits=c(-3.5, 3)) +
 # scale_y_continuous(limits=c(0, 1)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()






########## Not in Loops
ninloop <- read.table("HIVexp_NOTinLoopScored.txt", header=F, sep="\t")
colnames(ninloop) <- c("chrHIV","HIVstart","HIVend", "HIVbrcd", "HIVexp","strand", "chrom",	"start",	"end",	"H3K27ac",	"H3K36me3",	"H3K79me3",	"H3K9me3",	"H3K4me1",	"H3K4me3",	"RNAse",	"MNase",	"DNase",	"TTseq",	"H3K27me3")

### Separate by unique loops
unInLoop <- ninloop[!duplicated(ninloop[,c("chrHIV","HIVstart","HIVend")]),]
nunInLoop <- ninloop[duplicated(ninloop[,c("chrHIV","HIVstart","HIVend")]),]


### Plot unique by x = HIV expression y= MNase rpkm
pdf("BHIVE_expression_1x_inLoopscored.pdf")
ggplot(uInLoop, aes(x=HIVexp, y=MNase)) + 
  scale_color_manual(values=c("#0072B2")) +
  geom_point() + 
  #geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Host Gene Expression") + 
  xlab("Log10 HIV Expression") + 
  ylab("MNase-seq rpkm in loop") + 
  scale_x_continuous(limits=c(-3.5, 3)) +
  scale_y_continuous(limits=c(0, 1)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

### Plot unique by x = HIV expression y= MNase rpkm
pdf("BHIVE_expression_2plus_inLoopscored.pdf")
ggplot(nuInLoop, aes(x=HIVexp, y=MNase)) + 
  scale_color_manual(values=c("#0072B2")) +
  geom_point() + 
  #geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Host Gene Expression") + 
  xlab("Log10 HIV Expression") + 
  ylab("MNase-seq rpkm in loop") + 
  scale_x_continuous(limits=c(-3.5, 3)) +
  scale_y_continuous(limits=c(0, 1)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()









########## Tads
inloop <- read.table("HIVexp_inTadScored.txt", header=F, sep="\t")
colnames(inloop) <- c("chrHIV","HIVstart","HIVend", "HIVbrcd", "HIVexp","strand", "chrom",	"start",	"end",	"H3K27ac",	"H3K36me3",	"H3K79me3",	"H3K9me3",	"H3K4me1",	"H3K4me3",	"RNAse",	"MNase",	"DNase",	"TTseq",	"H3K27me3")

### Separate by unique loops
uInLoop <- inloop[!duplicated(inloop[,c("chrHIV","HIVstart","HIVend")]),]
nuInLoop <- inloop[duplicated(inloop[,c("chrHIV","HIVstart","HIVend")]),]

### Plot unique by x = HIV expression y= MNase rpkm
pdf("plots/BHIVE_expression_1x_inLoopscored.pdf")
ggplot(uInLoop, aes(x=HIVexp, y=H3K27me3)) + 
  scale_color_manual(values=c("#0072B2")) +
  geom_point() + 
  #geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Host Gene Expression") + 
  xlab("Log10 HIV Expression") + 
  ylab("MNase-seq rpkm in loop") + 
  scale_x_continuous(limits=c(-3.5, 3)) +
  #scale_y_continuous(limits=c(0, 1)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

### Plot unique by x = HIV expression y= MNase rpkm
pdf("plots/BHIVE_expression_2plus_inLoopscored.pdf")
ggplot(nuInLoop, aes(x=HIVexp, y=H3K27me3)) + 
  scale_color_manual(values=c("#0072B2")) +
  geom_point() + 
  #geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression to Host Gene Expression") + 
  xlab("Log10 HIV Expression") + 
  ylab("MNase-seq rpkm in loop") + 
  scale_x_continuous(limits=c(-3.5, 3)) +
 # scale_y_continuous(limits=c(0, 1)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()