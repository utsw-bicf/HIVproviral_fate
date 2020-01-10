### This R script makes a jitter scatter plot of HIV expression values for each of the 5 groups

library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(scales)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots')

########################################
### Load bed files by group
g1 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed", header =F, sep="\t")
g2 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed", header =F, sep="\t")
g3 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed", header =F, sep="\t")
g4 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed", header =F, sep="\t")
g5 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed", header =F, sep="\t")

### add groups
g1$group <- "Intergenic Same"
g2$group <- "Intergenic Convergent"
g3$group <- "Intergenic Divergent"
g4$group <- "Intragenic Same"
g5$group <- "Intragenic Convergent"

### Merge Tables
df <- rbind(g1, g2, g3, g4, g5)

### Make plot, jitter scatter plot
pdf("BHIVE_expression_by_group.pdf")
ggplot(df, aes(x=group, y=V11)) + 
  geom_jitter(position=position_jitter(width=.2, height=0)) +
  ggtitle("HIV Expression by Cluster") + 
  ylab("HIV Expression in log10") +
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.5, color="red") +
  scale_x_discrete(labels = wrap_format(10)) +
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

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots')

### Load file
g1 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed", header =F, sep="\t")
g2 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed", header =F, sep="\t")
g3 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed", header =F, sep="\t")
g4 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed", header =F, sep="\t")
g5 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed", header =F, sep="\t")

### add groups
g1$group <- "1"
g2$group <- "2"
g3$group <- "3"
g4$group <- "4"
g5$group <- "5"

### Merge Tables
df <- rbind(g1, g2, g3, g4, g5)
df2 <- df[c("V11", "group")]
colnames(df2) <- c("HIVexp", "group")

### Group 6
g6 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6thru8_2gene.bed", header =F, sep="\t")
g6$group <- "6"
df3 <- g6[c("V15","group")]
colnames(df3) <- c("HIVexp", "group")

### reduce to 2 important columsn
rEL <- rbind(df2,df3)

### make levels in proper order
rEL$group <- factor(rEL$group,levels = c("1", "2", "3", "4", "5", "6"))
levels(rEL$group)

### Calculate summary stats
library(dplyr)
sumStats <- group_by(rEL, group) %>%
  summarise(
    count = n(),
    mean = mean(HIVexp, na.rm = TRUE),
    sd = sd(HIVexp, na.rm = TRUE)
  )

# Compute the analysis of variance
res.aov <- aov(HIVexp ~ group, data = rEL)

# Summary of the analysis
summary(res.aov)

### This result means that there ARE differences in the group
### Need to do a two-way to determine which group

# Tukey multiple pairwise-comparisons
TukeyHSD(res.aov)

# Multiple comparisons using multcomp package
library(multcomp)
summary(glht(res.aov, linfct = mcp(group = "Tukey")))

# Pairewise t-test
pairwise.t.test(rEL$HIVexp, rEL$group, p.adjust.method = "BH")

# Check the homogeneity of variance assumption
plot(res.aov, 1)

library(car)
leveneTest(HIVexp ~ group, data = rEL)

### !!! This test was significant, cannot assume homogeneity of variances
# ANOVA test with no assumption of equal variances
oneway.test(HIVexp ~ group, data = rEL)

# Pairwise t-tests with no assumption of equal variances
pairwise.t.test(rEL$HIVexp, rEL$group,
                p.adjust.method = "BH", pool.sd = FALSE)

# Check the normality assumption
plot(res.aov, 2)

# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )

### !!! This test was significant, cannot assume normality 
kruskal.test(HIVexp ~ group, data = rEL)

pairwise.wilcox.test(rEL$HIVexp, rEL$group, p.adjust.method = "BH")
