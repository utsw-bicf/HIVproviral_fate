### This script makes a circos plot of HIV expression and sub insertion type

library(RCircos)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/circos')

### Load HIV exp data
group1 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed", header = F, sep = "\t")
group2 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed", header = F, sep = "\t")
group3 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed", header = F, sep = "\t")
group4 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed", header = F, sep = "\t")
group5 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed", header = F, sep = "\t")
group6 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6_2gene_opposite_noexp.bed", header = F, sep = "\t")
group7 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group7_2gene_same_noexp.bed", header = F, sep = "\t")
group8 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group8_2gene_opposite_noexp.bed", header = F, sep = "\t")

### Load reference ideogram
data("UCSC.HG38.Human.CytoBandIdeogram")

### Initialize RCircos Core Components
chr.exclude <- NULL
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
tracks.inside <- 1
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

### Plot circos
out.file <- "test.pdf"
pdf(file=out.file, height=8, width=8, compress=TRUE)
RCircos.Set.Plot.Area()
par(mai=c(0.25, 0.25, 0.25, 0.25))
plot.new()
plot.window(c(-2.5,2.5), c(-2.5, 2.5))

# get chromosomes
RCircos.Chromosome.Ideogram.Plot()

# add labels
name.col <- 4
side <- "in"
track.num <- 1


RCircos.Histogram.Plot(group8, data.col=4, track.num=4, side="in")

dev.off()


#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################

#install.packages("shiny")  
#install.packages("circlize")  
#install.packages("RColorBrewer")
#install.packages("data.table")
#install.packages("RLumShiny")  
#source("https://bioconductor.org/biocLite.R")  
#biocLite("GenomicRanges")

library("shiny")
library("circlize")
library("RColorBrewer")
library("data.table")
library("RLumShiny")
library("GenomicRanges")

run_as shiny

setwd('/work/BICF/s185797/programs/shinyCircos/shinyCircos')

# Load ui.R and server.R and run app

### Make dataframe
group1 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed", header = F, sep = "\t")
group2 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed", header = F, sep = "\t")
group3 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed", header = F, sep = "\t")
group4 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed", header = F, sep = "\t")
group5 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed", header = F, sep = "\t")
group6 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6_2gene_opposite_noexp.bed", header = F, sep = " ")
group7 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group7_2gene_same_noexp.bed", header = F, sep = " ")
group8 <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group8_2gene_opposite_noexp.bed", header = F, sep = " ")
group1$group <- 1
group2$group <- 2
group3$group <- 3
group4$group <- 4
group5$group <- 5
group6$group <- 6
group7$group <- 7
group8$group <- 8

group1 <- subset(group1, select = c(1,2,3,24))
group2 <- subset(group2, select = c(1,2,3,24))
group3 <- subset(group3, select = c(1,2,3,24))
group4 <- subset(group4, select = c(1,2,3,24))
group5 <- subset(group5, select = c(1,2,3,24))
group6 <- subset(group6, select = c(1,2,3,5))
group7 <- subset(group7, select = c(1,2,3,5))
group8 <- subset(group8, select = c(1,2,3,5))

df <- rbind(group1,group2,group3,group4,group5,group6,group7,group8)

write.table(df,file=("groupt.csv"), quote=F, col.names = F, row.names = F, sep=",")