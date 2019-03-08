### this R script makes the scatter plots for bhive results
library(ggplot2)
library(tidyr)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS')

################################################### Ivan's Scatter Plots ###################################################
### Load the same and opposite bed files
same <- read.table("hiv_expression_closestTSS_same_plus_dis.bed", header=F, sep="\t")
opp <- read.table("hiv_expression_closestTSS_opposite_plus_dis.bed", header=F, sep="\t")

### add log10(distance)
same$logDis <- log10(same$V23)
opp$logDis <- log10(opp$V23)

### Pearson's correlation
corsame <- cor.test(same$logDis, same$V11, method = "pearson", conf.level = 0.95)
coropp <- cor.test(opp$logDis, opp$V11, method = "pearson", conf.level = 0.95)


### Make scatter plot
### add pearson correlation manually
pdf("BHIVE_expression_scatterplot_closestTSS_same.pdf")
ggplot(same, aes(x=logDis, y=V11, size=V8)) + 
  geom_point(alpha=1/4) + 
  ggtitle("Same Orientation") + 
  xlab("Log10 Distance to TSS") + 
  ylab("HIV Expression") + 
  annotate("text", x = 6, y = 1, label = "r = -0.11") + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()


pdf("BHIVE_expression_scatterplot_closestTSS_opposite.pdf")
ggplot(opp, aes(x=logDis, y=V11, size=V8)) + 
  geom_point(alpha=1/4, color="red") + 
  ggtitle("Opposite Orientation") + 
  xlab("Log10 Distance to TSS") + 
  ylab("HIV Expression") + 
  annotate("text", x = 6, y = 1, label = "r = -0.13") + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()
