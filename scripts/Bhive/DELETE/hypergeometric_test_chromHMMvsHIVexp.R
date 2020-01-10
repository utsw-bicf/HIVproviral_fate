### This script does a hypergeometric test
### After looking at BHIVE_exp_by_chromstates_all.pdf, a majority of HIV insertions are in states 1,2,and 10
### Now we want to test if that insertions is by chance or statisticall over represented
### This is based on chromHMM state length vs total genome
library(IRanges)
library(broom)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_chromstates')


### Load BHIVE expression
hivexp <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt", header=T, sep='\t')

### Load chromHMM; 15 states, histone only
chromHMM <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed", header=F, sep='\t', skip=1)
chromHMM$V4 <- gsub("U", "\\1", chromHMM$V4)

### Add state to BHIVE expression
ihivexp <- with(hivexp, IRanges(locus, width=1, names=chr))
ichromHMMs <- with(chromHMM, IRanges(V2, V3, names=V1))
olaps <- findOverlaps(ihivexp, ichromHMMs)
df <- cbind(hivexp[queryHits(olaps),], chromHMM[subjectHits(olaps),])
df2 <- df[which(as.character(df$chr) == as.character(df$V1)),]


######################################################################
######################################################################
######################################################################
### Loop through hypergeomtric

total <- sum(chromHMM$V3-chromHMM$V2+1)
pr2 <- data.frame()
for (i in c(1:15)) {
  HMMr <- df2[df2$V4 == i,]
  chromr <- chromHMM[chromHMM$V4 == i,]
  HMMtotal <- sum(chromr$V3-chromr$V2+1)
  NROW(HMMr)
  NROW(df2)
  HMMtotal/total
  hgt <- as.data.frame(t(rbind(phyper(NROW(HMMr), HMMtotal, (total-HMMtotal), NROW(df2), lower.tail=F), dhyper(NROW(HMMr), HMMtotal, (total-HMMtotal), NROW(df2), log=T))))
  hgt$state <- i
  hgt$actual <- NROW(HMMr)
  hgt$percentgenome <- format(HMMtotal/total,5)
  hgt$perHMM <- (NROW(HMMr))/(NROW(df2))
  pr2 <- rbind(pr2,hgt)
}
#rownames(pr2) <- c(1:15)
colnames(pr2) <- c("phyper", "log(dhyper)", "state", "insertion_count", "state_pct_genome", "pct_insertions")

write.table(pr2, file=("hypergeometric_test.txt"), quote=F, row.names = F, sep='\t')

########## Loop through states, chance
total <- sum(chromHMM$V3-chromHMM$V2+1)
pr1 <- data.frame()
for (i in c(1:15)) {
  HMMr <- df2[df2$V4 == i,]
  chromr <- chromHMM[chromHMM$V4 == i,]
  HMMtotal <- sum(chromr$V3-chromr$V2+1)
  NROW(HMMr)
  NROW(df2)
  HMMtotal/total
  pt <- tidy(prop.test(x=c((NROW(HMMr)),HMMtotal),n=c((NROW(df2)),total),p=NULL,alternative="two.sided", conf.level=0.95, correct=F))
  pt$state <- i
  pt$actual <- NROW(HMMr)
  pt$percentgenome <- format(HMMtotal/total,5)
  pt$perHMM <- (NROW(HMMr))/(NROW(df2))
  pr1 <- rbind(pr1,pt)
}

##########
########## Loop through states, more
total <- sum(chromHMM$V3-chromHMM$V2+1)
prm <- data.frame()
for (i in c(1:15)) {
  HMMr <- df2[df2$V4 == i,]
  chromr <- chromHMM[chromHMM$V4 == i,]
  HMMtotal <- sum(chromr$V3-chromr$V2+1)
  NROW(HMMr)
  NROW(df2)
  HMMtotal/total
  pt <- tidy(prop.test(x=c((NROW(HMMr)),HMMtotal),n=c((NROW(df2)),total),p=NULL,alternative="greater", conf.level=0.95, correct=F))
  pt$state <- i
  pt$actual <- NROW(HMMr)
  pt$percentgenome <- format(HMMtotal/total,5)
  pt$perHMM <- (NROW(HMMr))/(NROW(df2))
  prm <- rbind(prm,pt)
}

##########
########## Loop through states, less
total <- sum(chromHMM$V3-chromHMM$V2+1)
prl <- data.frame()
for (i in c(1:15)) {
  HMMr <- df2[df2$V4 == i,]
  chromr <- chromHMM[chromHMM$V4 == i,]
  HMMtotal <- sum(chromr$V3-chromr$V2+1)
  NROW(HMMr)
  NROW(df2)
  HMMtotal/total
  pt <- tidy(prop.test(x=c((NROW(HMMr)),HMMtotal),n=c((NROW(df2)),total),p=NULL,alternative="less", conf.level=0.95, correct=F))
  pt$state <- i
  pt$actual <- NROW(HMMr)
  pt$percentgenome <- format(HMMtotal/total,5)
  pt$perHMM <- (NROW(HMMr))/(NROW(df2))
  prl <- rbind(prl,pt)
}

write.table(pr1, file=("prob.test_two.ways.txt"), quote=F, row.names = F, sep='\t')
write.table(prl, file=("prob.test_less.txt"), quote=F, row.names = F, sep='\t')
write.table(prm, file=("prob.test_more.txt"), quote=F, row.names = F, sep='\t')


