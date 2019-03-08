### This script identifies insertions where the closest TSS is not the closest Gene

setwd("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated")

### Load the two files
TSS <- read.table("hiv_expression_intergenic_closestTSS.bed", header=F, sep="\t")
Gene <- read.table("hiv_expression_intergenic_closestGene.bed", header=F, sep="\t")

### merge the two files
df <- merge(TSS, Gene, by=(c("V1", "V2", "V3")), all=T)

same <- df[as.character(df$V12.x)==as.character(df$V12.y),]
dif <- df[as.character(df$V12.x)!=as.character(df$V12.y),]

### Write merge tables
write.table(same,file="hiv_expression_intergenic_closestTSS_eq_closestGene.bed", quote=F,row.names=F,col.names=F, sep='\t')
write.table(dif,file="hiv_expression_intergenic_closestTSS_noteq_closestGene.bed", quote=F,row.names=F,col.names=F, sep='\t')