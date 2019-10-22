### This script makes a pie chart and
### jitter plot of HIV expression vs repeat masker
library(ggplot2)
library(RColorBrewer)
library(scales)


setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/HIVexpression_vs_repeats')

hr <- read.table("HIV_expression_repeats_fixed.bed", header = F, sep = "\t")

### Add group number
g1 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same_latency.bed", header = F, sep = "\t")
g2 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite_latency.bed", header = F, sep = "\t")
g3 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite_latency.bed", header = F, sep = "\t")
g4 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene_latency.bed", header = F, sep = "\t")
g5 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene_latency.bed", header = F, sep = "\t")
g6 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6thru8_2gene.bed", header = F, sep = "\t")

g1r <- subset(g1,select = c(V1,V2))
g1r$group <- as.character("group1")
colnames(g1r) <- c("chr", "brcd", "group")

g2r <- subset(g2,select = c(V1,V2))
g2r$group <- as.character("group2")
colnames(g2r) <- c("chr", "brcd", "group")

g3r <- subset(g3,select = c(V1,V2))
g3r$group <- as.character("group3")
colnames(g3r) <- c("chr", "brcd", "group")

g4r <- subset(g4,select = c(V1,V2))
g4r$group <- as.character("group4")
colnames(g4r) <- c("chr", "brcd", "group")

g5r <- subset(g5,select = c(V1,V2))
g5r$group <- as.character("group5")
colnames(g5r) <- c("chr", "brcd", "group")

g6r <- subset(g6,select = c(V1,V8))
g6r$group <- as.character("group6")
colnames(g6r) <- c("chr", "brcd", "group")

g <- rbind(g1r, g2r, g3r, g4r, g5r, g6r)

hrg <- merge(hr, g, by.x=c("V1", "V4"), by.y=c("chr", "brcd"))

pdf("HIVexp_by_repeatmasker.pdf")
ggplot(hrg, aes(x=V10, y=V5)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5, aes(color = group)) +
  ggtitle("HIV Expression by RepeatMasker") + 
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
        legend.position="right")
dev.off()

##### Make pie chart
pdf("piechart_HIVexp_by_repeatmasker.pdf")
mytable <- table(hrg$V10)
lbls <- paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls,
    main="HIV location by RepeatMasker")
dev.off