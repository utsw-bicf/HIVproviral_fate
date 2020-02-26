########## This file recreates Bhive vs Lamin Subcompartments, but with Jurkat
########## Try A vs B, then try A1, A2, B
library(IRanges)
library(broom)
library(ggplot2)
library(scales)

########## A1, A2 vs B
setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_lamin/jurkat')

### Load BHIVE expression
hivexp <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt", header=T, sep='\t')

### Load Lamin; jurkat A1,A2,B
A <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1_kmeans.bed", header=F, sep='\t')
B <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1.bed", header=F, sep='\t')
Lamin <- rbind(A,B)
Lamin$L <- gsub("Lamin","", Lamin$V4)

### Add state to BHIVE expression
ihivexp <- with(hivexp, IRanges(locus, width=1, names=chr))
iLamins <- with(Lamin, IRanges(V2, V3, names=V1))
olaps <- findOverlaps(ihivexp, iLamins)
df <- cbind(hivexp[queryHits(olaps),], Lamin[subjectHits(olaps),])
df2 <- df[which(as.character(df$chr) == as.character(df$V1)),]

# Remove the ones with NA's
df2 <- df2[!df2$brcd == "CTCTTTTCCTCGGAA", ]

### Plot
pdf("BHIVE_exp_by_lamin_A1A2vsB_jurkat.pdf")
ggplot(df2, aes(x=L, y=expr)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
  ggtitle("HIV Expression by Lamin subcompartments") + 
  ylab("HIV Expression in log10") +
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.5, color="red") +
  scale_x_discrete(labels = wrap_format(10)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        #axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.x = element_text(hjust = 1),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

### Do the same but with a jitter box-plot
pdf("BHIVE_exp_by_lamin_A1A2vsB_jurkat_boxplot.pdf")
ggplot(df2, aes(x=L, y=expr, fill=L)) + 
  geom_boxplot() +
  #geom_jitter(position=position_jitter(width=.3, height=0), size= 0.2) +
  ggtitle("HIV Expression by Lamin subcompartments") + 
  ylab("HIV Expression in log10") +
  scale_fill_manual(values=c("#FF0000","#006400", "#007FFF")) +
  scale_x_discrete(labels = wrap_format(10)) +
  scale_y_continuous(limits=c(-4, 4)) +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()


#### Hypergeometric test
total <- sum(Lamin$V3-Lamin$V2+1)
pr2 <- data.frame()
for (i in c("LaminA1","LaminA2","LaminB")) {
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

colnames(pr2) <- c("phyper", "log(dhyper)", "state", "insertion_count", "state_pct_genome", "pct_insertions")

write.table(pr2, file=("Lamin_A1A2vsB_hypergeometric_test_jurkat.txt"), quote=F, row.names = F, sep='\t')


### Kruskal-Wallis rank sum test
### make levels in proper order
rEL <- df2[,c("expr", "L")]

rEL$L <- factor(rEL$L,levels = c("A1", "A2", "B"))
levels(rEL$L)

### Calculate summary stats
library(dplyr)
sumStats <- group_by(rEL, L) %>%
  summarise(
    count = n(),
    mean = mean(expr, na.rm = TRUE),
    sd = sd(expr, na.rm = TRUE)
  )

# Compute the analysis of variance
res.aov <- aov(expr ~ L, data = rEL)

# Summary of the analysis
summary(res.aov)

### This result means that there ARE differences in the group
### Need to do a two-way to determine which group

# Tukey multiple pairwise-comparisons
TukeyHSD(res.aov)

# Multiple comparisons using multcomp package
library(multcomp)
summary(glht(res.aov, linfct = mcp(L = "Tukey")))

# Pairewise t-test
pairwise.t.test(rEL$expr, rEL$L, p.adjust.method = "BH")

# Check the homogeneity of variance assumption
plot(res.aov, 1)

library(car)
leveneTest(expr ~ L, data = rEL)

### !!! This test was not significant, we can assume homogeneity of variances
## ANOVA test with no assumption of equal variances
#oneway.test(expr ~ L, data = rEL)

## Pairwise t-tests with no assumption of equal variances
#pairwise.t.test(rEL$expr, rEL$L,
#                p.adjust.method = "BH", pool.sd = FALSE)

# Check the normality assumption
plot(res.aov, 2)

# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )

### !!! This test was significant, cannot assume normality 
kruskal.test(expr ~ L, data = rEL)

pairwise.wilcox.test(rEL$expr, rEL$L, p.adjust.method = "BH")





###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
########## A vs B
#setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_lamin/jurkat')

### Load BHIVE expression
#hivexp <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt", header=T, sep='\t')

### Load Lamin; 5 states plus NA
#Lamin <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/Jurkat_PC1_Lamin.bed", header=F, sep='\t')


### Add state to BHIVE expression
#ihivexp <- with(hivexp, IRanges(locus, width=1, names=chr))
#iLamins <- with(Lamin, IRanges(V2, V3, names=V1))
#olaps <- findOverlaps(ihivexp, iLamins)
#df <- cbind(hivexp[queryHits(olaps),], Lamin[subjectHits(olaps),])
#df2 <- df[which(as.character(df$chr) == as.character(df$V1)),]

# Remove the ones with NA's
#df2 <- df2[!df2$brcd == "CTCTTTTCCTCGGAA", ]

### Plot
#pdf("BHIVE_exp_by_lamin_AvsB_jurkat.pdf")
#ggplot(df2, aes(x=V4, y=expr)) + 
#  geom_violin(trim = FALSE) +
#  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
#  ggtitle("HIV Expression by Lamin subcompartments") + 
#  ylab("HIV Expression in log10") +
#  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.5, color="red") +
#  scale_x_discrete(labels = wrap_format(10)) +
#  theme(plot.title = element_text(size=22, hjust = 0.5), 
#        #axis.text.x = element_text(angle = 90, hjust = 1),
#        axis.text.x = element_text(hjust = 1),
#        axis.title.x=element_blank(),
#        panel.grid.major = element_blank(), 
#        panel.grid.minor = element_blank(),
#        panel.background = element_blank(), 
#        axis.line = element_line(colour = "black"),
#        legend.position="none")
#dev.off()

#### Hypergeometric test
#total <- sum(Lamin$V3-Lamin$V2+1)
#pr2 <- data.frame()
#for (i in c("LaminA","LaminB")) {
#  HMMr <- df2[df2$V4 == i,]
#  laminr <- Lamin[which(Lamin$V4 == i),]
#  HMMtotal <- sum(laminr$V3-laminr$V2+1)
#  NROW(HMMr)
#  NROW(df2)
#  HMMtotal/total
#  hgt <- as.data.frame(t(rbind(phyper(NROW(HMMr), HMMtotal, (total-HMMtotal), NROW(df2), lower.tail=F), dhyper(NROW(HMMr), HMMtotal, (total-HMMtotal), NROW(df2), log=T))))
#  hgt$state <- i
#  hgt$actual <- NROW(HMMr)
#  hgt$percentgenome <- format(HMMtotal/total,5)
#  hgt$perHMM <- (NROW(HMMr))/(NROW(df2))
#  pr2 <- rbind(pr2,hgt)
#}

#colnames(pr2) <- c("phyper", "log(dhyper)", "state", "insertion_count", "state_pct_genome", "pct_insertions")

#write.table(pr2, file=("Lamin_AvsB_hypergeometric_test_jurkat.txt"), quote=F, row.names = F, sep='\t')


###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
########## A1, A2, B1, B2, B3
#setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_lamin/jurkat')

### Load BHIVE expression
#hivexp <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt", header=T, sep='\t')

### Load Lamin; 5 states plus NA
#Lamin <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/Jurkat_PC1withKmeans_Lamin_A1A2B1B2B3.bed", header=F, sep='\t')


### Add state to BHIVE expression
#ihivexp <- with(hivexp, IRanges(locus, width=1, names=chr))
#iLamins <- with(Lamin, IRanges(V2, V3, names=V1))
#olaps <- findOverlaps(ihivexp, iLamins)
#df <- cbind(hivexp[queryHits(olaps),], Lamin[subjectHits(olaps),])
#df2 <- df[which(as.character(df$chr) == as.character(df$V1)),]

# Remove the ones with NA's
#df2 <- df2[!df2$brcd == "CTCTTTTCCTCGGAA", ]

### Plot
#pdf("BHIVE_exp_by_lamin_A1A2B1B2B3_jurkat.pdf")
#ggplot(df2, aes(x=V4, y=expr)) + 
#  geom_violin(trim = FALSE) +
#  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
#  ggtitle("HIV Expression by Lamin subcompartments") + 
#  ylab("HIV Expression in log10") +
#  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.5, color="red") +
#  scale_x_discrete(labels = wrap_format(10)) +
#  theme(plot.title = element_text(size=22, hjust = 0.5), 
#        #axis.text.x = element_text(angle = 90, hjust = 1),
#        axis.text.x = element_text(hjust = 1),
#        axis.title.x=element_blank(),
#        panel.grid.major = element_blank(), 
#        panel.grid.minor = element_blank(),
#        panel.background = element_blank(), 
#        axis.line = element_line(colour = "black"),
#        legend.position="none")
#dev.off()

#### Hypergeometric test
#total <- sum(Lamin$V3-Lamin$V2+1)
#pr2 <- data.frame()
#for (i in c("LaminA1","LaminA2","LaminB1","LaminB2","LaminB3")) {
#  HMMr <- df2[df2$V4 == i,]
#  laminr <- Lamin[which(Lamin$V4 == i),]
#  HMMtotal <- sum(laminr$V3-laminr$V2+1)
#  NROW(HMMr)
#  NROW(df2)
#  HMMtotal/total
#  hgt <- as.data.frame(t(rbind(phyper(NROW(HMMr), HMMtotal, (total-HMMtotal), NROW(df2), lower.tail=F), dhyper(NROW(HMMr), HMMtotal, (total-HMMtotal), NROW(df2), log=T))))
#  hgt$state <- i
#  hgt$actual <- NROW(HMMr)
#  hgt$percentgenome <- format(HMMtotal/total,5)
#  hgt$perHMM <- (NROW(HMMr))/(NROW(df2))
#  pr2 <- rbind(pr2,hgt)
#}

#colnames(pr2) <- c("phyper", "log(dhyper)", "state", "insertion_count", "state_pct_genome", "pct_insertions")

#write.table(pr2, file=("Lamin_A1A2B1B2B3_hypergeometric_test_jurkat.txt"), quote=F, row.names = F, sep='\t')

### Kruskal-Wallis rank sum test
### make levels in proper order
#rEL <- df2[,c("expr", "V4")]

#rEL$V4 <- factor(rEL$V4,levels = c("LaminA1", "LaminA2", "LaminB1", "LaminB2", "LaminB3"))
#levels(rEL$V4)

### Calculate summary stats
#library(dplyr)
#sumStats <- group_by(rEL, V4) %>%
#  summarise(
#    count = n(),
#    mean = mean(expr, na.rm = TRUE),
#    sd = sd(expr, na.rm = TRUE)
#  )

# Compute the analysis of variance
#res.aov <- aov(expr ~ V4, data = rEL)

# Summary of the analysis
#summary(res.aov)

### This result means that there ARE differences in the group
### Need to do a two-way to determine which group

# Tukey multiple pairwise-comparisons
#TukeyHSD(res.aov)

# Multiple comparisons using multcomp package
#library(multcomp)
#summary(glht(res.aov, linfct = mcp(V4 = "Tukey")))

# Pairewise t-test
#pairwise.t.test(rEL$expr, rEL$V4, p.adjust.method = "BH")

# Check the homogeneity of variance assumption
#plot(res.aov, 1)

#library(car)
#leveneTest(expr ~ V4, data = rEL)

### !!! This test was significant, cannot assume homogeneity of variances
# ANOVA test with no assumption of equal variances
#oneway.test(expr ~ V4, data = rEL)

# Pairwise t-tests with no assumption of equal variances
#pairwise.t.test(rEL$expr, rEL$V4,
#                p.adjust.method = "BH", pool.sd = FALSE)

# Check the normality assumption
#plot(res.aov, 2)

# Extract the residuals
#aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test
#shapiro.test(x = aov_residuals )

### !!! This test was significant, cannot assume normality 
#kruskal.test(expr ~ V4, data = rEL)

#pairwise.wilcox.test(rEL$expr, rEL$V4, p.adjust.method = "BH")


