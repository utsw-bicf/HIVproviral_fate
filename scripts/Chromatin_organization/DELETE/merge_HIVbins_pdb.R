##### This script adds HIV bins to pdb bins

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000')

pdb <- read.table("MDS.HiC_combine_20000.pdb", header = F)
colnames(pdb) <- c("atom", "bin", "value", "edg", "a", "value1", "value2", "value3", "value4", "value5")
hivb_high <- read.table("HIV_high_bin.bed", header = F)
hivb_med <- read.table("HIV_med_bin.bed", header = F)
hivb_low <- read.table("HIV_low_bin.bed", header = F)

### Merge together and output
pdb_high <- merge(hivb_high, pdb, by.x = "V10", by.y = "bin", all.x = T)
pdb_med <- merge(hivb_med, pdb, by.x = "V10", by.y = "bin", all.x = T)
pdb_low <- merge(hivb_low, pdb, by.x = "V10", by.y = "bin", all.x = T)

dist <- function(x){
  sqrt((x - 100 )^ 2 + (y - 100 )^ 2 + (z - 100 )^ 2)
}

pdb_high$dist <- sqrt((pdb_high$value1-100)^2 + (pdb_high$value2-100)^2 + (pdb_high$value3-100)^2)
pdb_med$dist <- sqrt((pdb_med$value1-100)^2 + (pdb_med$value2-100)^2 + (pdb_med$value3-100)^2)
pdb_low$dist <- sqrt((pdb_low$value1-100)^2 + (pdb_low$value2-100)^2 + (pdb_low$value3-100)^2)

write.table(pdb_high, file=("HIV_high_bin_pdb.bed"), quote = F, sep = "\t", col.names = F, row.names = F)
write.table(pdb_med, file=("HIV_med_bin_pdb.bed"), quote = F, sep = "\t", col.names = F, row.names = F)
write.table(pdb_low, file=("HIV_low_bin_pdb.bed"), quote = F, sep = "\t", col.names = F, row.names = F)

min(pdb_high$dist)
max(pdb_high$dist)
mean(pdb_high$dist)
median(pdb_high$dist)

min(pdb_med$dist)
max(pdb_med$dist)
mean(pdb_med$dist)
median(pdb_med$dist)


min(pdb_low$dist)
max(pdb_low$dist)
mean(pdb_low$dist)
median(pdb_low$dist)

