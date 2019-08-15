### This script makes a pie chart that represents the size of the lamin subcompartments

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_lamin')
Lamin <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments_hg38.bed", header=F, sep='\t')

A1 <- Lamin[which(Lamin$V4 == "A1"),]
A1total <- sum(A1$V3-A1$V2+1)

A2 <- Lamin[which(Lamin$V4 == "A2"),]
A2total <- sum(A2$V3-A2$V2+1)

B1 <- Lamin[which(Lamin$V4 == "B1"),]
B1total <- sum(B1$V3-B1$V2+1)

B2 <- Lamin[which(Lamin$V4 == "B2"),]
B2total <- sum(B2$V3-B2$V2+1)

B3 <- Lamin[which(Lamin$V4 == "B3"),]
B3total <- sum(B3$V3-B3$V2+1)

B4 <- Lamin[which(Lamin$V4 == "B4"),]
B4total <- sum(B4$V3-B4$V2+1)

Nann <- 3088286401-A1total-A2total-B1total-B2total-B3total-B4total


slices <- c(A1total,A2total,B1total,B2total,B3total,B4total,Nann) 
lbls <- c("A1","A2","B1","B2","B3","B4","NA")
pdf("Lamin_subcompartments_pie.pdf")
pie(slices, labels=lbls, main="Lamin Subcompartments")
dev.off()

total <- Nann+A1total+A2total+B1total+B2total+B3total+B4total
A1total/total
A2total/total
B1total/total
B2total/total
B3total/total
B4total/total
Nann/total

