### This script adds expression values to rpkm output

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning')

# Load HIV
HIV <- read.table("input_files/hiv_expression_table.txt", header=F, sep="\t")
colnames(HIV) <- c("chrom","start","end", "HIVexp")

# Load both rpkm outputs
r1 <- read.table("HIVexpression_histones_filtered_peaks.tsv", header=T, sep="\t")
r2 <- read.table("HIVexpression_HisPlus3_filtered_peaks.tsv", header=T, sep="\t")

# add expression
r1e <- merge(HIV, r1, by.x=c("chrom","start","end"), by.y=c("chrom","start","end"), all.y = T)
r2e <- merge(HIV, r2, by.x=c("chrom","start","end"), by.y=c("chrom","start","end"), all.y = T)

# Write out
write.table(r1e, file="HIVexpression_histones_4ML.tsv", quote=F, col.names=T, sep="\t")
write.table(r2e, file="HIVexpression_HisPlus3_4ML.tsv", quote=F, col.names=T, sep="\t")