### This R script makes a violin scatter plot by x=Lamin y=HIVexp
### Then does a hypergeometric test

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_lamin')

### Load file
expLamin <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/older_test_files/HIVexp_HisPlus4_chrom_lamin_4ML.tsv", header = T, sep = "\t")

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
  #ggtitle("HIV Expression by Lamin subcompartments") + 
  ylab("HIV Expression in log10") +
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.5, color="red") +
  scale_x_discrete(labels = wrap_format(10)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        #axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

########## pie chart of inserts
expLamin <- expLamin[!expLamin$chr == "chrUn_GL000195v1", ]

### reduce to 2 important columsn
rEL <- expLamin[c("HIVexp", "lamin_0_200")]
rEL$lamin_0_200[is.na(rEL$lamin_0_200)] <- "Unknown"

### make levels in proper order
rEL$lamin_0_200 <- factor(rEL$lamin_0_200,levels = c("A1", "A2", "B1", "B2", "B3", "B4", "Unknown"))

pdf("piechart_HIVinserts_GSMlamin.pdf")
mytable <- table (rEL$lamin_0_200)
lbls <- paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls, main="Insertions in Sub-compartments")
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


######################################################################
######################################################################
######################################################################
######################################################################
#### One way annova
#### followed: http://www.sthda.com/english/wiki/one-way-anova-test-in-r

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_lamin')

### Load file
expLamin <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/older_test_files/HIVexp_His7Plus4_chrom_lamin_4ML.tsv", 
                       header = T, sep = "\t", stringsAsFactors = FALSE)

expLamin <- expLamin[!expLamin$chr == "chrUn_GL000195v1", ]

### reduce to 2 important columsn
rEL <- expLamin[c("HIVexp", "lamin_0_200")]
rEL$lamin_0_200[rEL$lamin_0_200 == '0'] <- "Unknown"

### make levels in proper order
rEL$lamin_0_200 <- factor(rEL$lamin_0_200,levels = c("A1", "A2", "B1", "B2", "B3", "B4", "Unknown"))
levels(rEL$lamin_0_200)

### Calculate summary stats
library(dplyr)
sumStats <- group_by(rEL, lamin_0_200) %>%
  summarise(
    count = n(),
    mean = mean(HIVexp, na.rm = TRUE),
    sd = sd(HIVexp, na.rm = TRUE)
  )

# Compute the analysis of variance
res.aov <- aov(HIVexp ~ lamin_0_200, data = rEL)

# Summary of the analysis
summary(res.aov)

### This result means that there ARE differences in the group
### Need to do a two-way to determine which group

# Tukey multiple pairwise-comparisons
TukeyHSD(res.aov)

# Multiple comparisons using multcomp package
library(multcomp)
summary(glht(res.aov, linfct = mcp(lamin_0_200 = "Tukey")))

# Pairewise t-test
pairwise.t.test(rEL$HIVexp, rEL$lamin_0_200, p.adjust.method = "BH")

# Check the homogeneity of variance assumption
plot(res.aov, 1)

library(car)
leveneTest(HIVexp ~ lamin_0_200, data = rEL)

### !!! This test was significant, cannot assume homogeneity of variances
# ANOVA test with no assumption of equal variances
oneway.test(HIVexp ~ lamin_0_200, data = rEL)

# Pairwise t-tests with no assumption of equal variances
pairwise.t.test(rEL$HIVexp, rEL$lamin_0_200,
                p.adjust.method = "BH", pool.sd = FALSE)

# Check the normality assumption
plot(res.aov, 2)

# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )

### !!! This test was significant, cannot assume normality 
kruskal.test(HIVexp ~ lamin_0_200, data = rEL)

pairwise.wilcox.test(rEL$HIVexp, rEL$lamin_0_200, p.adjust.method = "BH")
