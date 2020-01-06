########## Runs Kmeans clustering of Lamin subcompartments

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged')

LA <- read.table("LaminA_positive_PC1.bed", header = F, sep = "\t")

hist(LA$V5)

set.seed(20)
test <- kmeans(LA[,5], 2)
LA$K <- test$cluster

write.table(LA, file=("LaminA_positive_PC1_kmeans.txt"), quote = F, col.names = F, row.names = F, sep = "\t")


LB <- read.table("LaminB_negative_PC1.bed", header = F, sep = "\t")

hist(LB$V5)

set.seed(20)
test2 <- kmeans(LB[,5], 3)
LB$K <- test2$cluster

write.table(LB, file=("LaminB_negative_PC1_kmeans_3.txt"), quote = F, col.names = F, row.names = F, sep = "\t")


test3 <- kmeans(LB[,5], 4)
LB$K <- test3$cluster

write.table(LB, file=("LaminB_negative_PC1_kmeans_4.txt"), quote = F, col.names = F, row.names = F, sep = "\t")



###### Lamin B, inparticular B4 on chr19 is different, so
###### Kmeans cluster all but chr19 with 3 groups
###### Kmeans cluster chr19 with 4 groups
LB <- read.table("LaminB_negative_PC1.bed", header = F, sep = "\t")

LBnot19 <- LB[which(LB$V1 != "chr19"),]
LB19 <- LB[which(LB$V1 == "chr19"),]

hist(LBnot19$V5)
hist(LB19$V5)

set.seed(20)
test10 <- kmeans(LBnot19[,5], 3)
LBnot19$K <- test10$cluster

set.seed(20)
test11 <- kmeans(LB19[,5],4)
LB19$K <- test11$cluster

write.table(LBnot19, file=("LaminB_negative_PC1_kmeans3_notChr19.txt"), quote = F, col.names = F, row.names = F, sep = "\t")
write.table(LB19, file=("LaminB_negative_PC1_kmeans4_Chr19.txt"), quote = F, col.names = F, row.names = F, sep = "\t")


### chr19 only has 3 clusters, try kmeans 3
LB <- read.table("LaminB_negative_PC1.bed", header = F, sep = "\t")
LB19 <- LB[which(LB$V1 == "chr19"),]
set.seed(20)
test11 <- kmeans(LB19[,5],3)
LB19$K <- test11$cluster
write.table(LB19, file=("LaminB_negative_PC1_kmeans3_Chr19.txt"), quote = F, col.names = F, row.names = F, sep = "\t")
