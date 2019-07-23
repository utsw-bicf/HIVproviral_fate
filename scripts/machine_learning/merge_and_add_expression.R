### This script merges the 2 kb by 200 bp files
### chromHMM, HisPlus4, lamin, HIVexpression

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files')
### Load files and merge by barcode
HIVexp <- read.table("hiv_expression_table_withname.txt", header=F, sep="\t")
df <- HIVexp
colnames(df) <- c("chr", "start", "end", "bc", "HIVexp")

### Add chromHMM
chromHMM <- Sys.glob("hiv_expression_ML_*_chromHMM.bed")
for (i in 1:length(chromHMM)) {
  rd <- read.csv(file = chromHMM[i], header = F, stringsAsFactors = F, sep="\t")
  rd <- rd[,c(1,4,5)]
  fn <- sub(pattern = "hiv_expression_ML_(.*)_chromHMM.bed", replacement = "\\1", basename(chromHMM[i]))
  colnames(rd) <- c("chr","bc", paste("chromHMM", fn, sep="_"))
  df <- merge(df, rd, by = c("bc","chr"), all.x = T)
}

### Add lamin
lamin <- Sys.glob("hiv_expression_ML_*_lamin.bed")
for (i in 1:length(lamin)) {
  rd2 <- read.csv(file = lamin[i], header = F, stringsAsFactors = F, sep="\t")
  rd2 <- rd2[,c(1,4,5)]
  fn2 <- sub(pattern = "hiv_expression_ML_(.*)_lamin.bed", replacement = "\\1", basename(lamin[i]))
  colnames(rd2) <- c("chr","bc", paste("lamin", fn2, sep="_"))
  df <- merge(df, rd2, by = c("bc","chr"), all.x = T)
}
df[is.na(df)] <- 0

### Add Histone + 4 data
HisP4 <- Sys.glob("hiv_expression_ML_*HisPlus4_filtered_peaks.tsv")
for (i in 1:length(HisP4)) {
  rd3 <- read.csv(file = HisP4[i], header = T, stringsAsFactors = F, sep="\t")
  rd3 <- rd3[,-c(2:3)]
  fn3 <- sub(pattern = "hiv_expression_ML_(.*)_HisPlus4_filtered_peaks.tsv", replacement = "\\1", basename(HisP4[i]))
  colnames(rd3) <- paste(colnames(rd3),fn3, sep = "_")
  colnames(rd3)[1] <- "chr"
  colnames(rd3)[2] <- "bc"
  df <- merge(df, rd3, by = c("chr","bc"), all.x = T)
}

### Make missing data 0
df[is.na(df)] <- 0

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/')
write.table(df, file=("HIVexp_HisPlus4_chrom_lamin_4ML.tsv"), quote=F, col.names = T, row.names = F, sep = "\t")