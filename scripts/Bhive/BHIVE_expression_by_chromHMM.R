### This R script makes a violin plot
### X = chromHMM state
### Y = HIV expression

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_chromstates')


### Load BHIVE expression
hivexp <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt", header=T, sep='\t')

### Load chromHMM; 15 states, 7 histones histone
chromHMM <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed", header=F, sep='\t', skip=1)

### Add state to BHIVE expression
library(IRanges)
ihivexp <- with(hivexp, IRanges(locus, width=1, names=chr))
ichromHMMs <- with(chromHMM, IRanges(V2, V3, names=V1))
olaps <- findOverlaps(ihivexp, ichromHMMs)
df <- cbind(hivexp[queryHits(olaps),], chromHMM[subjectHits(olaps),])
df2 <- df[which(as.character(df$chr) == as.character(df$V1)),]

### remove chr from column chr, adding to new column
df2$cnum <- gsub("chr", "\\1", df2$chr)
### remove U from chromHMM state
df2$HMM <- gsub("U", "\\1", df2$V4)

### Plot
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(scales)

df2$HMM <- factor(df2$HMM,levels = c("1", "2","3","4","5","6","7","8","9","10", "11","12","13","14","15"))

pdf("BHIVE_exp_by_chromstates_histone7_Ernstorder.pdf")
ggplot(df2, aes(x=HMM, y=expr)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0)) +
  ggtitle("chromHMM state of HIV expression") + 
  ylab("HIV Expression") +
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.5, color="red") +
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



######################################################################
######################################################################
######################################################################
######################################################################
#### One way annova
#### followed: http://www.sthda.com/english/wiki/one-way-anova-test-in-r

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_chromstates')

### Load file
hivexp <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt", header=T, sep='\t')

### Load chromHMM; 15 states, 7 histones histone
chromHMM <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed", header=F, sep='\t', skip=1)

### Add state to BHIVE expression
library(IRanges)
ihivexp <- with(hivexp, IRanges(locus, width=1, names=chr))
ichromHMMs <- with(chromHMM, IRanges(V2, V3, names=V1))
olaps <- findOverlaps(ihivexp, ichromHMMs)
df <- cbind(hivexp[queryHits(olaps),], chromHMM[subjectHits(olaps),])
df2 <- df[which(as.character(df$chr) == as.character(df$V1)),]

### remove chr from column chr, adding to new column
df2$cnum <- gsub("chr", "\\1", df2$chr)
### remove U from chromHMM state
df2$HMM <- gsub("U", "\\1", df2$V4)

df2 <- df2[!df2$chr == "chrUn_GL000195v1", ]

### reduce to 2 important columsn
rEL <- df2[c("expr", "HMM")]

### make levels in proper order
rEL$HMM <- factor(rEL$HMM,levels = c("1", "2","3","4","5","6","7","8","9","10", "11","12","13","14","15"))
levels(rEL$HMM)

### Calculate summary stats
library(dplyr)
sumStats <- group_by(rEL, HMM) %>%
  summarise(
    count = n(),
    mean = mean(expr, na.rm = TRUE),
    sd = sd(expr, na.rm = TRUE)
  )

# Compute the analysis of variance
res.aov <- aov(expr ~ HMM, data = rEL)

# Summary of the analysis
summary(res.aov)

### This result means that there ARE differences in the group
### Need to do a two-way to determine which group

# Tukey multiple pairwise-comparisons
TukeyHSD(res.aov)

# Multiple comparisons using multcomp package
library(multcomp)
summary(glht(res.aov, linfct = mcp(HMM = "Tukey")))

# Pairewise t-test
pairwise.t.test(rEL$expr, rEL$HMM, p.adjust.method = "BH")

# Check the homogeneity of variance assumption
plot(res.aov, 1)

library(car)
leveneTest(expr ~ HMM, data = rEL)

### !!! This test was significant, cannot assume homogeneity of variances
# ANOVA test with no assumption of equal variances
oneway.test(expr ~ HMM, data = rEL)

# Pairwise t-tests with no assumption of equal variances
pairwise.t.test(rEL$expr, rEL$HMM,
                p.adjust.method = "BH", pool.sd = FALSE)

# Check the normality assumption
plot(res.aov, 2)

# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )

### !!! This test was significant, cannot assume normality 
kruskal.test(expr ~ HMM, data = rEL)

pairwise.wilcox.test(rEL$expr, rEL$HMM, p.adjust.method = "BH")
