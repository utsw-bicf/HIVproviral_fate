### This R script makes a violin plot of
### X = latent 
### Y = nearest RNA expression

setwd("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles")

### Load RNAseq data
fpkm <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq_reg/workflow/output_PE/countTable.fpkm.txt", header=T, sep="\t")
fpkm$log10fpkm <- log(rowMeans(fpkm[c('Emily_jurkat_A', 'Emily_jurkat_B', 'Emily_jurkat_C')], na.rm=TRUE)+1)

### Load the file with latencies; 1: insert is latent (GFP-), 0: insert is GFP active (GFP+)
latent <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/BHIVE_latent_edited.txt", header = T, sep ="\t", skip=21)
### reduce latent to chr, barcode, latent
latentR <- latent[,c("chrom", "brcd", "latent")]
latentR$latent[latentR$latent == 0] <- "GFP+"
latentR$latent[latentR$latent == 1] <- "GFP-"
latentR$latent[latentR$latent == 2] <- "both"

### Load the annotation files; this identifies the nearest TSS to each insertion
group1 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed", header = F, sep = "\t")
group2 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed", header = F, sep = "\t")
group3 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed", header = F, sep = "\t")
group4 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed", header = F, sep = "\t")
group5 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed", header = F, sep = "\t")

catgroup <- rbind(group1,group2,group3,group4,group5)

### Make a new dataframe with latent and gene identities
latentLoc <- merge(catgroup, latentR, by.x=c("V1","V4"), by.y=c("chrom","brcd"), all=F)
### Reduce dataframe and make correct ensembl number
latentLoc$ENSEMBL <- substr(latentLoc$V15, 0, 15)
latentLocR <- latentLoc[,c("ENSEMBL", "latent")]

### Merge RNAseq and latent
expLatenat <- merge(latentLocR, fpkm, by="ENSEMBL", all=F)

### Plot
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(scales)

expLatenat$latent <- factor(expLatenat$latent,levels = c("GFP-","GFP+","both"))
pdf("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_latency/Latency_by_RNAseqExp.pdf")
ggplot(expLatenat, aes(x=latent, y=log10fpkm)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0)) +
  ggtitle("Latency") + 
  ylab("Host Gene Expression (Log10(fpkm))") +
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
