### This script adds whether the expression is latent based on the barcode and chromosome.
### Got the latent file from the publication.

setwd("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles")
### Load the annotation files
group1 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed", header = F, sep = "\t")
group2 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed", header = F, sep = "\t")
group3 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed", header = F, sep = "\t")
group4 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed", header = F, sep = "\t")
group5 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed", header = F, sep = "\t")

### Load the file with latencies; 1: insert is latent (GFP-), 0: insert is GFP active (GFP+)
latent <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/BHIVE_latent_edited.txt", header = T, sep ="\t", skip=21)

### reduce latent to chr, barcode, latent
latentR <- latent[,c("chrom", "brcd", "latent")]
latentR$latent[latentR$latent == 0] <- "GFP+"
latentR$latent[latentR$latent == 1] <- "GFP-"
latentR$latent[latentR$latent == 2] <- "both"

### merge tables together
group1m <- merge(group1, latentR, by.x=c("V1","V4"), by.y=c("chrom","brcd"), all.x=T)
group2m <- merge(group2, latentR, by.x=c("V1","V4"), by.y=c("chrom","brcd"), all.x=T)
group3m <- merge(group3, latentR, by.x=c("V1","V4"), by.y=c("chrom","brcd"), all.x=T)
group4m <- merge(group4, latentR, by.x=c("V1","V4"), by.y=c("chrom","brcd"), all.x=T)
group5m <- merge(group5, latentR, by.x=c("V1","V4"), by.y=c("chrom","brcd"), all.x=T)

### Make "NA" in latent to unknown
group1m$latent[is.na(group1m$latent)] <- "unknown"
group2m$latent[is.na(group2m$latent)] <- "unknown"
group3m$latent[is.na(group3m$latent)] <- "unknown"
group4m$latent[is.na(group4m$latent)] <- "unknown"
group5m$latent[is.na(group5m$latent)] <- "unknown"

### Write out
write.table(group1m, file=("group1_Intergenic_same_latency.bed"), quote=F, col.names = F, row.names = F, sep='\t')
write.table(group2m, file=("group2_Intergenic_downstream_opposite_latency.bed"), quote=F, col.names = F, row.names = F, sep='\t')
write.table(group3m, file=("group3_Intergenic_upstream_opposite_latency.bed"), quote=F, col.names = F, row.names = F, sep='\t')
write.table(group4m, file=("group4_Intragenic_same_1gene_latency.bed"), quote=F, col.names = F, row.names = F, sep='\t')
write.table(group5m, file=("group5_Intragenic_opposite_1gene_latency.bed"), quote=F, col.names = F, row.names = F, sep='\t')


######################################################################
######################################################################
######################################################################
########## Chart latency (by group) in x vs. expression in y; violin plot
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(scales)

group1m$latent <- factor(group1m$latent,levels = c("GFP-","GFP+","both","unknown"))
pdf("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_latency/group1_Intergenic_same_latency.pdf")
ggplot(group1m, aes(x=latent, y=V11)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0)) +
  ggtitle("Group 1 latency") + 
  ylab("HIV Expression in log10") +
  xlab("GFP Status") +
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


group2m$latent <- factor(group2m$latent,levels = c("GFP-","GFP+","both","unknown"))
pdf("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_latency/group2_Intergenic_downstream_opposite_latency.pdf")
ggplot(group2m, aes(x=latent, y=V11)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0)) +
  ggtitle("Group 2 latency") + 
  ylab("HIV Expression in log10") +
  xlab("GFP Status") +
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


group3m$latent <- factor(group3m$latent,levels = c("GFP-","GFP+","both","unknown"))
pdf("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_latency/group3_Intergenic_upstream_opposite_latency.pdf")
ggplot(group3m, aes(x=latent, y=V11)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0)) +
  ggtitle("Group 3 latency") + 
  ylab("HIV Expression in log10") +
  xlab("GFP Status") +
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


group4m$latent <- factor(group4m$latent,levels = c("GFP-","GFP+","both","unknown"))
pdf("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_latency/group4_Intragenic_same_1gene_latency.pdf")
ggplot(group4m, aes(x=latent, y=V11)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0)) +
  ggtitle("Group 4 latency") + 
  ylab("HIV Expression in log10") +
  xlab("GFP Status") +
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


group5m$latent <- factor(group5m$latent,levels = c("GFP-","GFP+","both","unknown"))
pdf("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_latency/group5_Intragenic_opposite_1gene_latency.pdf")
ggplot(group5m, aes(x=latent, y=V11)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0)) +
  ggtitle("Group 5 latency") + 
  ylab("HIV Expression in log10") +
  xlab("GFP Status") +
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