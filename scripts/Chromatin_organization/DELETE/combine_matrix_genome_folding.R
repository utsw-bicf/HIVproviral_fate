########## This script combines the matrix files from Hi-C pro output
########## use raw

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix')

dir.create("combine")
dir.create("combine/raw")

sizes=c('10000','20000','40000','150000','500000','1000000')
#x='1000000'
for (x in sizes) {
  f1 <- read.table(file=(paste0("rep1/raw/",x,"/rep1_",x,".matrix")), header=F,sep="\t")
  f2 <- read.table(file=(paste0("rep2/raw/",x,"/rep2_",x,".matrix")), header=F,sep="\t")
  c <- merge(f1, f2, by = c("V1","V2"), all = T)
  c$count <- c$V3.x +c$V3.y
  nc <- c[,c("V1","V2","count")]
  
  dir.create(paste0("combine/raw/",x))
  write.table(nc, file=paste("combine/raw/",x,"/combine_",x,".matrix", sep=""), quote = F, sep = "\t", col.names = F, row.names = F)
}
