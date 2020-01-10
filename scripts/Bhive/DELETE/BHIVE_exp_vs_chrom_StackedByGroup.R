### This R script makes a violin scatter plot by x=chromosome, y=Hiv expression
### But divided by group (1-6)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots')

library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(scales)
#library(cowplot)
library(grid)
require(gridExtra)

### Load data by group
group1 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed", header = F, sep = "\t")
group2 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed", header = F, sep = "\t")
group3 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed", header = F, sep = "\t")
group4 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed", header = F, sep = "\t")
group5 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed", header = F, sep = "\t")
#group6 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6_2gene_opposite_noexp.bed", header = F, sep = " ")
#group7 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group7_2gene_same_noexp.bed", header = F, sep = "\t")
#group8 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group8_2gene_opposite_noexp.bed", header = F, sep = "\t")
#group6to8 <- rbind(group6,group7,group8)

### Get rest
hivexp <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt", header=T, sep='\t')
g1to5 <- rbind(group1,group2,group3, group4, group5)
gb <- g1to5$V4
hb <- hivexp$brcd
wantlist <- hb[which(!(hb %in% gb))]
group6to8 <- filter(hivexp, hivexp$brcd %in% wantlist)
group6to8 <- group6to8[!(group6to8$chr == "chrUn_GL000195v1"),]

### Set chromosomes
group1$V1 <- factor(group1$V1,levels = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                           "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                                           "chr21","chr22","chrX","chrY"))
group2$V1 <- factor(group2$V1,levels = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                         "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                                         "chr21","chr22","chrX","chrY"))
group3$V1 <- factor(group3$V1,levels = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                         "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                                         "chr21","chr22","chrX","chrY"))
group4$V1 <- factor(group4$V1,levels = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                         "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                                         "chr21","chr22","chrX","chrY"))
group5$V1 <- factor(group5$V1,levels = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                         "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                                         "chr21","chr22","chrX","chrY"))
group6to8$chr <- factor(group6to8$chr,levels = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                         "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                                         "chr21","chr22","chrX","chrY"))


### Make idividual plots
pdf("BHIVE_exp_by_chrom_StackedByGroup.pdf")
g1p <- ggplot(group1, aes(x=V1, y=V11)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
  #ggtitle("HIV Expression by Chromosome") + 
  ylab("Group 1") +
  scale_x_discrete(labels = wrap_format(10), drop=FALSE) +
  scale_y_continuous(limits=c(-4.5,3)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        #axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")

g2p <- ggplot(group2, aes(x=V1, y=V11)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
  #ggtitle("HIV Expression by Chromosome") + 
  ylab("Group 2") +
  scale_x_discrete(labels = wrap_format(10), drop=FALSE) +
  scale_y_continuous(limits=c(-4.5,3)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")

g3p <- ggplot(group3, aes(x=V1, y=V11)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
  #ggtitle("HIV Expression by Chromosome") + 
   ylab("Group 3") +
  scale_x_discrete(labels = wrap_format(10), drop=FALSE) +
  scale_y_continuous(limits=c(-4.5,3)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")

g4p <- ggplot(group4, aes(x=V1, y=V11)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
  #ggtitle("HIV Expression by Chromosome") + 
   ylab("Group 4") +
  scale_x_discrete(labels = wrap_format(10), drop=FALSE) +
  scale_y_continuous(limits=c(-4.5,3)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")

g5p <- ggplot(group5, aes(x=V1, y=V11)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
  #ggtitle("HIV Expression by Chromosome") + 
   ylab("Group 5") +
  scale_x_discrete(labels = wrap_format(10), drop=FALSE) +
  scale_y_continuous(limits=c(-4.5,3)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")

g6p <- ggplot(group6to8, aes(x=chr, y=expr)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
  #ggtitle("HIV Expression by Chromosome") + 
   ylab("Group 6") +
  scale_x_discrete(labels = wrap_format(10), drop=FALSE) +
  scale_y_continuous(limits=c(-4.5,3)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")

#title <- ggdraw() + draw_label("HIV Expression by Chromosome", fontface='bold')

#grid.newpage()
#grid.draw(rbind(ggplotGrob(title),ggplotGrob(g1p), ggplotGrob(g2p), 
#                ggplotGrob(g3p), ggplotGrob(g4p), ggplotGrob(g5p),
#                ggplotGrob(g6p),size = "last"))
#                #ggplotGrob(g6p),size = c(1,2,2,2,2,2)))
grid.arrange(g1p, g2p, g3p, g4p, g5p, g6p, ncol = 1, 
             heights = c(2,2,2,2,2,3),
             top = textGrob("HIV Expression by Chromosome", vjust = 1, 
                            gp = gpar(fontface = "bold", cex = 1.5)))
#grid.arrange(g1p, g2p, g3p, g4p, g5p, g6p, ncol = 2, heights = c(2,2,2,2,2,2),
#             top = textGrob("HIV Expression by Chromosom", vjust = 1, 
#                            gp = gpar(fontface = "bold", cex = 1.5)))

dev.off()



###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################


setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots')

library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(scales)
library(grid)
require(gridExtra)

### Load data by group
group1 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed", header = F, sep = "\t")
group2 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed", header = F, sep = "\t")
group3 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed", header = F, sep = "\t")
group4 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed", header = F, sep = "\t")
group5 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed", header = F, sep = "\t")

### Get rest
hivexp <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt", header=T, sep='\t')
g1to5 <- rbind(group1,group2,group3, group4, group5)
gb <- g1to5$V4
hb <- hivexp$brcd
wantlist <- hb[which(!(hb %in% gb))]
group6to8 <- filter(hivexp, hivexp$brcd %in% wantlist)
group6to8 <- group6to8[!(group6to8$chr == "chrUn_GL000195v1"),]

### Rename chromosomes
group1$V1 <- gsub("chr", "\\1", group1$V1)
group2$V1 <- gsub("chr", "\\1", group2$V1)
group3$V1 <- gsub("chr", "\\1", group3$V1)
group4$V1 <- gsub("chr", "\\1", group4$V1)
group5$V1 <- gsub("chr", "\\1", group5$V1)
group6to8$chr <- gsub("chr", "\\1", group6to8$chr)

### Set chromosomes
group1$V1 <- factor(group1$V1,levels = c("1", "2","3","4","5","6","7","8","9","10",
                                         "11","12","13","14","15","16","17","18","19","20",
                                         "21","22","X","Y"))
group2$V1 <- factor(group2$V1,levels = c("1", "2","3","4","5","6","7","8","9","10",
                                         "11","12","13","14","15","16","17","18","19","20",
                                         "21","22","X","Y"))
group3$V1 <- factor(group3$V1,levels = c("1", "2","3","4","5","6","7","8","9","10",
                                         "11","12","13","14","15","16","17","18","19","20",
                                         "21","22","X","Y"))
group4$V1 <- factor(group4$V1,levels = c("1", "2","3","4","5","6","7","8","9","10",
                                         "11","12","13","14","15","16","17","18","19","20",
                                         "21","22","X","Y"))
group5$V1 <- factor(group5$V1,levels = c("1", "2","3","4","5","6","7","8","9","10",
                                         "11","12","13","14","15","16","17","18","19","20",
                                         "21","22","X","Y"))
group6to8$chr <- factor(group6to8$chr,levels = c("1", "2","3","4","5","6","7","8","9","10",
                                                 "11","12","13","14","15","16","17","18","19","20",
                                                 "21","22","X","Y"))


### Make idividual plots
pdf("BHIVE_exp_by_chrom_StackedByGroup2.pdf")
g1p <- ggplot(group1, aes(x=V1, y=V11)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
  #ggtitle("HIV Expression by Chromosome") + 
  ylab("log10 HIV Expression") +
 # xlab("Chromosome Number") +
  scale_x_discrete(drop=FALSE) +
  scale_y_continuous(limits=c(-4.5,4.5)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        #axis.text.x = element_text(angle = 90, hjust = 1),
        #axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")

g2p <- ggplot(group2, aes(x=V1, y=V11)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
  #ggtitle("HIV Expression by Chromosome") + 
 # ylab("log10 HIV Expression") +
 # xlab("Chromosome Number") +
  scale_x_discrete(labels = wrap_format(10), drop=FALSE) +
  scale_y_continuous(limits=c(-4.5,4.5)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")

g3p <- ggplot(group3, aes(x=V1, y=V11)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
  #ggtitle("HIV Expression by Chromosome") + 
 # ylab("log10 HIV Expression") +
 # xlab("Chromosome Number") +
  scale_x_discrete(labels = wrap_format(10), drop=FALSE) +
  scale_y_continuous(limits=c(-4.5,4.5)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")

g4p <- ggplot(group4, aes(x=V1, y=V11)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
  #ggtitle("HIV Expression by Chromosome") + 
  #ylab("log10 HIV Expression") +
  #xlab("Chromosome Number") +
  scale_x_discrete(labels = wrap_format(10), drop=FALSE) +
  scale_y_continuous(limits=c(-4.5,4.5)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")

g5p <- ggplot(group5, aes(x=V1, y=V11)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
  #ggtitle("HIV Expression by Chromosome") + 
 # ylab("log10 HIV Expression") +
  #xlab("Chromosome Number") +
  scale_x_discrete(labels = wrap_format(10), drop=FALSE) +
  scale_y_continuous(limits=c(-4.5,4.5)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")

g6p <- ggplot(group6to8, aes(x=chr, y=expr)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
  #ggtitle("HIV Expression by Chromosome") + 
 # ylab("log10 HIV Expression") +
  #xlab("Chromosome Number") +
  scale_x_discrete(labels = wrap_format(10), drop=FALSE) +
  scale_y_continuous(limits=c(-4.5,4.5)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        #axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")

grid.arrange(g1p, g2p, g3p, g4p, g5p, g6p, ncol = 1, 
             heights = c(2,2,2,2,2,2))


dev.off()