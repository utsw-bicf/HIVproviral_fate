### This script makes the list
### Use in script BHIVE_corr_2insertions_1gene.sh

setwd("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots")


g4 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed", header =F, sep="\t")
g5 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed", header =F, sep="\t")

### Merge Tables
dfintra <- rbind(g4, g5)

### Find multiples in 1 gene
n_occur <- data.frame(table(dfintra$V15))

### Work only with 2 insertions in 1 gene; >2 will be handled differently
list <- n_occur[n_occur$Freq == 2,]
list2 <- dfintra[dfintra$V15 %in% list$Var1,]
slist <- list2[order(list2$V15),]
#slist$rep <- c(1,2)
write.table(slist, file=("delete_sorted_list.txt"), row.names = F, col.names = F, sep="\t", quote =F)



