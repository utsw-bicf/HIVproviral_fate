### This script makes a pie chart that represents the size of the lamin subcompartments

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_lamin/jurkat')

### Load Lamin; jurkat A1,A2,B
A <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1_kmeans.bed", header=F, sep='\t')
B <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1.bed", header=F, sep='\t')
Lamin <- rbind(A,B)
Lamin$L <- gsub("Lamin","", Lamin$V4)


A1 <- Lamin[which(Lamin$L == "A1"),]
A1total <- sum(A1$V3-A1$V2+1)

A2 <- Lamin[which(Lamin$L == "A2"),]
A2total <- sum(A2$V3-A2$V2+1)

B1 <- Lamin[which(Lamin$L == "B"),]
B1total <- sum(B1$V3-B1$V2+1)

Nann <- 3088286401-A1total-A2total-B1total


slices <- c(A1total,A2total,B1total,Nann) 
lbls <- c("A1","A2","B", "NA")
pdf("Lamin_subcompartments_jurkat_pie.pdf")
pie(slices, labels=lbls, main="Lamin Subcompartments")
dev.off()

total <- Nann+A1total+A2total+B1total
A1total/total
A2total/total
B1total/total
Nann/total

