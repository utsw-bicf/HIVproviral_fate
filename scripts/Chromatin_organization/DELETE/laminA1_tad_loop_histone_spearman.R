### This script takes the input of 
### 1) Lamin A1 only
### 2) A loop in Lamin A1
### 3) A Tad in Lamin A1
### calculates the spearman correlation of histon marks (RPKM of 1-3 area) for each HIV insertion
### Sorted by HIV

library("ggpubr")
library("ggcorrplot")

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/')

### Load in RPKM
A1 <- read.table("Lamin_filtered_peaks.tsv", header = T, sep = "\t")
loop <- read.table("Loop_filtered_peaks.tsv", header = T, sep = "\t")
tad <- read.table("Tad_filtered_peaks.tsv", header = T, sep = "\t")

### Add HIV expression
HIV <- read.table("HIVexp_A1lamin_inloopuniq_intad.bed", header = F, sep = "\t")
HIV <- HIV[,c(1:5)]
colnames(HIV) <- c("chrom", "HIVstart", "HIVend", "name", "exp")

### Merge tables and sort high to low HIV exp
A1m <- merge(A1,HIV, by = c("chrom", "name"), all = T)
A1m <- A1m[order(-(A1m$exp)),]
loopm <- merge(loop,HIV, by = c("chrom", "name"), all = T)
loopm <- loopm[order(-(loopm$exp)),]
tadm <- merge(tad,HIV, by = c("chrom", "name"), all = T)
tadm <- tadm[order(-(tadm$exp)),]


########### What kind of test can we run?
### Plot values
ggscatter(A1m, x = "exp", y = "H3K27ac", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "HIV expression", ylab = "RPKM H3K27ac in Lamin A1/Loop/Tad")

### Is the data normal distribution
### Result: if pvalue > 0.05 we can assume normality
A1mexpST <- shapiro.test(A1m$exp)
A1mH3K27acST <- shapiro.test(A1m$H3K27ac)
A1mH3K36me3ST <- shapiro.test(A1m$H3K36me3)
A1mH3K79me3ST <- shapiro.test(A1m$H3K79me3)
A1mH3K9me3ST <- shapiro.test(A1m$H3K9me3)
A1mH3K4me1ST <- shapiro.test(A1m$H3K4me1)
A1mH3K4me3ST <- shapiro.test(A1m$H3K4me3)
A1mRNAseST <- shapiro.test(A1m$RNAse)
A1mMNaseST <- shapiro.test(A1m$MNase)
A1mDNaseST <- shapiro.test(A1m$DNase)
A1mTTseqST <- shapiro.test(A1m$TTseq)

ggpubr::ggqqplot(A1m$exp, ylab = "HIV expression")
ggpubr::ggqqplot(A1m$H3K27ac, ylab = "RPKM H3K27ac in Lamin A1/Loop/Tad")
# Note: this data is not normally distributed based on pvalue < 0.05
# This means that I cannot use a Pearson correlation, but will have to use a spearman or kendell




########## Kendall and spearman
A1mK <- cor.test(A1m$exp, A1m$H3K27ac,  method="kendall")
A1mS <- cor.test(A1m$exp, A1m$H3K27ac,  method="spearman")

# If pvalue < 0.05, result: strong case for not correlated

########## Loop through all tests
LR <- data.frame(matrix(vector(),0,3))
### Lamin
for (i in 5:15) {
  a <- cor.test(A1m$exp, A1m[,i], method="kendall")
  print(paste(colnames(A1m)[i], " est:", a$estimate, " p=value:", a$p.value))
  #nd <- data.frame(paste(colnames(A1m)[i], a$estimate, a$p.value))
  #LR <- rbind(LR,nd)
  LR <- rbind(LR,data.frame(colnames(A1m)[i], a$estimate, a$p.value))
  
  pdf(file=paste0("kendall_corr_plots/", colnames(A1m)[i], "_LaminA1.pdf"))
  p <- ggscatter(A1m, x = "exp", y = colnames(A1m)[i], 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "kendall",
                 add.params = list(color = "blue", fill = "lightgray"),
                 xlab = "HIV expression", ylab = "RPKM in Lamin A1 containing HIV in Loop and Tad")
  print(p)
  dev.off()
}



### Make n.s. as NA
LR$a.estimate[LR$a.p.value > 0.05] <- "NA"
### Make the reverse of pvalue
LR$size <- 0.05-LR$a.p.value
LR$size[LR$a.p.value > 0.05] <- 0
LR$size[LR$a.p.value <= 0.05] <- 1
LR$cor <- abs(LR$a.estimate)

ggballoonplot(LR, x = 1, y = "colnames.A1m..i.", fill = "a.estimate", size = "size") +
  scale_fill_viridis_c(option = "C")

ggballoonplot(LR, x = 1, y = "colnames.A1m..i.", fill = "size", size = "cor", size.range = c(0, 1))

# [1] "H3K27ac  est: 0.123197809550348  p=value: 0.0208371613621535"
# [1] "H3K36me3  est: -0.0697477254068744  p=value: 0.190766451600759"
# [1] "H3K79me3  est: 0.136328610568254  p=value: 0.0105510275815696"
# [1] "H3K9me3  est: -0.0975541275624389  p=value: 0.0672645148169768"
# [1] "H3K4me1  est: -0.0473481236704474  p=value: 0.374461431999895"
# [1] "H3K4me3  est: 0.120571649346767  p=value: 0.0237183556113591"
# [1] "RNAse  est: 0.148378051502332  p=value: 0.00538175808034558"
# [1] "MNase  est: 0.0708290854907019  p=value: 0.183980280048569"
# [1] "DNase  est: 0.134093866938424  p=value: 0.0118960459400614"
# [1] "TTseq  est: 0.0949279673588578  p=value: 0.0749710359746672"
# [1] "H3K27me3  est: -0.0726828456344062  p=value: 0.172765439015291"
### 

### Loop
for (i in 5:15) {
  a <- cor.test(loopm$exp, loopm[,i], method="kendall")
  print(paste(colnames(loopm)[i], " est:", a$estimate, " p=value:", a$p.value))
  
  pdf(file=paste0("kendall_corr_plots/", colnames(loopm)[i], "_LaminA1_Loop.pdf"))
  p <- ggscatter(loopm, x = "exp", y = colnames(loopm)[i], 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "kendall",
                 add.params = list(color = "blue", fill = "lightgray"),
                 xlab = "HIV expression", ylab = "RPKM in Loop containing HIV in Lamin A1 and Tad")
  print(p)
  dev.off()
}

### Tad
for (i in 5:15) {
  a <- cor.test(tadm$exp, tadm[,i], method="kendall")
  print(paste(colnames(tadm)[i], " est:", a$estimate, " p=value:", a$p.value))
  
  pdf(file=paste0("kendall_corr_plots/", colnames(tadm)[i], "_LaminA1_Loop_Tad.pdf"))
  p <- ggscatter(tadm, x = "exp", y = colnames(tadm)[i], 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "kendall",
                 add.params = list(color = "blue", fill = "lightgray"),
                 xlab = "HIV expression", ylab = "RPKM in Tad containing HIV in Lamin A1 and Loop")
  print(p)
  dev.off()
}