### This R script makes a violin scatter plot by x=Lamin y=HIVexp
### Then does a hypergeometric test

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_lamin')

### Load file
expLamin <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/HIVexp_HisPlus4_chrom_lamin_4ML.tsv", header = T, sep = "\t")

### Remove unknown chr
expLamin <- expLamin[!expLamin$chr == "chrUn_GL000195v1", ]

### Make plot
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(scales)

expLamin$lamin_0_200 <- factor(expLamin$lamin_0_200,levels = c("A1", "A2", "B1", "B2", "B3", "B4", "unknown"))

pdf("BHIVE_exp_by_lamin_subcompartments.pdf")
ggplot(expLamin, aes(x=lamin_0_200, y=HIVexp)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
  ggtitle("HIV Expression by Lamin subcompartments") + 
  ylab("HIV Expression in log10") +
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.5, color="red") +
  scale_x_discrete(labels = wrap_format(10)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

######################################################################
######################################################################
######################################################################
######################################################################
#### Hypergeometric test
library(IRanges)
library(broom)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_lamin')

### Load BHIVE expression
hivexp <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt", header=T, sep='\t')

### Load Lamin; 5 states plus NA
Lamin <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments_hg38.bed", header=F, sep='\t')

### Add state to BHIVE expression
ihivexp <- with(hivexp, IRanges(locus, width=1, names=chr))
iLamins <- with(Lamin, IRanges(V2, V3, names=V1))
olaps <- findOverlaps(ihivexp, iLamins)
df <- cbind(hivexp[queryHits(olaps),], Lamin[subjectHits(olaps),])
df2 <- df[which(as.character(df$chr) == as.character(df$V1)),]

# Remove the ones with NA's
df2 <- df2[!df2$brcd == "CTCTTTTCCTCGGAA", ]

######################################################################

########## Loop through states, chance
total <- sum(Lamin$V3-Lamin$V2+1)
pr1 <- data.frame()
for (i in c("A1","A2","B1","B2","B3","B4")) {
  HMMr <- df2[df2$V4 == i,]
  laminr <- Lamin[which(Lamin$V4 == i),]
  HMMtotal <- sum(laminr$V3-laminr$V2+1)
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
rownames(pr1) <- c("A1","A2","B1","B2","B3","B4")



##########
########## Loop through states, more
total <- sum(Lamin$V3-Lamin$V2+1)
prm <- data.frame()
for (i in c("A1","A2","B1","B2","B3","B4")) {
  HMMr <- df2[df2$V4 == i,]
  laminr <- Lamin[which(Lamin$V4 == i),]
  HMMtotal <- sum(laminr$V3-laminr$V2+1)
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
rownames(prm) <- c("A1","A2","B1","B2","B3","B4")

##########
########## Loop through states, less
total <- sum(Lamin$V3-Lamin$V2+1)
prl <- data.frame()
for (i in c("A1","A2","B1","B2","B3","B4")) {
  HMMr <- df2[df2$V4 == i,]
  laminr <- Lamin[which(Lamin$V4 == i),]
  HMMtotal <- sum(laminr$V3-laminr$V2+1)
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
rownames(prl) <- c("A1","A2","B1","B2","B3","B4")



######################################################################
### Loop through hypergeomtric

total <- sum(Lamin$V3-Lamin$V2+1)
pr2 <- data.frame()
for (i in c("A1","A2","B1","B2","B3","B4")) {
  HMMr <- df2[df2$V4 == i,]
  laminr <- Lamin[which(Lamin$V4 == i),]
  HMMtotal <- sum(laminr$V3-laminr$V2+1)
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
rownames(pr2) <- c("A1","A2","B1","B2","B3","B4")
colnames(pr2) <- c("phyper", "log(dhyper)", "state", "insertion_count", "state_pct_genome", "pct_insertions")


write.table(pr1, file=("Lamin_prob.test_two.ways.txt"), quote=F, row.names = F, sep='\t')
write.table(prl, file=("Lamin_prob.test_less.txt"), quote=F, row.names = F, sep='\t')
write.table(prm, file=("Lamin_prob.test_more.txt"), quote=F, row.names = F, sep='\t')
write.table(pr2, file=("Lamin_hypergeometric_test.txt"), quote=F, row.names = F, sep='\t')