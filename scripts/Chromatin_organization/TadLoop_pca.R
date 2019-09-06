library(ggplot2)
library(ggbiplot)
library(pheatmap)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_TadLoop')

########## Loops
inloop <- read.table("HIVexp_inLoopScored.txt", header=F, sep="\t")
colnames(inloop) <- c("chrHIV","HIVstart","HIVend", "HIVbrcd", "HIVexp","strand", "chrom",	"start",	"end",	"H3K27ac",	"H3K36me3",	"H3K79me3",	"H3K9me3",	"H3K4me1",	"H3K4me3",	"RNAse",	"MNase",	"DNase",	"TTseq",	"H3K27me3")

### Make label for HIV expression
inloop$amount <- ifelse(abs(inloop$HIVexp) < 0.3, "nc", 
                    ifelse(inloop$HIVexp >= 0.3,"high", "low"))

### Separate by unique loops
uInLoop <- inloop[!duplicated(inloop[,c("chrHIV","HIVstart","HIVend")]),]
nuInLoop <- inloop[duplicated(inloop[,c("chrHIV","HIVstart","HIVend")]),]

test2 <- uInLoop$amount


test <- prcomp(uInLoop[,10:20], center = T, scale. = T)
plot(test, type= "l")

plot2 <- ggbiplot(test, obs.scale = 1, var.scale = 1, 
                  groups = test2, ellipse = TRUE, 
                  circle = TRUE)
print(plot2)

test <- dist(t(uInLoop[,10:20]), method = 'euclidean')
test2 <- hclust(test, method = "complete")
plot(test2)
pheatmap(test2)

mat <- as.matrix(uInLoop[,10:20])
STREE <- hclust(dist(t(mat)))
zscores <- scale(t(STREE))

t1 <- data.matrix(t(uInLoop[,10:20]))
STREE <- hclust(dist(t(t1)))
zscores <- scale(t(STREE))
pheatmap(t1)


d1 <- scale(t(uInLoop[,10:20]))
d2 <- dist(d1, method = "euclidean")
d3 <- hclust(d2, method = "complete")




a1 <- scale(uInLoop[,10:20])
pheatmap(a1)


test <- dist(uInLoop[,10:20], method = 'euclidean')
test2 <- hclust(test, method = "complete")
plot(test2)


