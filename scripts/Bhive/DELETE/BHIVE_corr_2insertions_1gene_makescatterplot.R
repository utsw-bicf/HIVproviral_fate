### This script makes the scatter plot
### Use in script BHIVE_corr_2insertions_1gene.sh
library(ggplot2)

setwd("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots")

#####################################################################################################
########## All dots together
### Load in file
ex_val <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots/two_insertions_1gene.txt", header=F, sep="\t")
ex_val$ensembl <- substr(ex_val$V1, 0, 15)

### Add gene names
gn <- read.table("/project/shared/bicf_workflow_ref/human/GRCh38/genenames.txt", header=T, sep="\t")
gn_r <- gn[,c("ensembl","symbol")]

ex_val_gn <- merge(ex_val, gn_r, by="ensembl")

### Pearson correlation
cor_ex_val_gn <- cor.test(ex_val_gn$V2, ex_val_gn$V3, method = "pearson", conf.level = 0.95)


### Make scatterplot
pdf("BHIVE_expression_2insertions_1gene.pdf")
ggplot(ex_val_gn, aes(x=V2, y=V3)) + 
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression of 2 Insertions into 1 Gene") + 
  xlab("Insertion 1 (closest to TSS) Expression") + 
  ylab("Insertion 2 (closest to TES) Expression") + 
  annotate("text", x = -2.5, y = 1.1, label = "r = -0.04 (p-value = 0.67)", size=3, hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()


#####################################################################################################
########## 3 plots same, opposite, mixed compared to gene direction
same <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots/two_insertions_1gene_same.txt", header=F, sep="\t")
opp <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots/two_insertions_1gene_opposite.txt", header=F, sep="\t")
mix <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_closestTSS_scatterplots/two_insertions_1gene_mixed.txt", header=F, sep="\t")

same$ensembl <- substr(same$V1, 0, 15)
opp$ensembl <- substr(opp$V1, 0, 15)
mix$ensembl <- substr(mix$V1, 0, 15)

### Add gene names
gn <- read.table("/project/shared/bicf_workflow_ref/human/GRCh38/genenames.txt", header=T, sep="\t")
gn_r <- gn[,c("ensembl","symbol")]

same_gn <- merge(same, gn_r, by="ensembl")
opp_gn <- merge(opp, gn_r, by="ensembl")
mix_gn <- merge(mix, gn_r, by="ensembl")

### add groups and concatenate tables
same_gn$group <- "Same"
opp_gn$group <- "Opposite"
mix_gn$group <- "Mixed"

df <- rbind(same_gn, opp_gn, mix_gn)

### Pearson correlation
cor_same_gn <- cor.test(same_gn$V2, same_gn$V3, method = "pearson", conf.level = 0.95)
cor_opp_gn <- cor.test(opp_gn$V2, opp_gn$V3, method = "pearson", conf.level = 0.95)
cor_mix_gn <- cor.test(mix_gn$V2, mix_gn$V3, method = "pearson", conf.level = 0.95)

### Make scatter plot by group
### add pearson correlation manually
pdf("BHIVE_expression_2insertions_1gene_bydirection.pdf")
ggplot(df, aes(x=V2, y=V3, color=group)) + 
  scale_color_manual(values=c("#CC79A7", "#D55E00", "#0072B2")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression of 2 Insertions into 1 Gene") + 
  xlab("Insertion 1 (closest to TSS) Expression") + 
  ylab("Insertion 2 (closest to TES) Expression") + 
  scale_x_continuous(limits=c(-2.5, 2.5)) +
  scale_y_continuous(limits=c(-2.5, 2)) +
  annotate("text", x = -2.5, y = 2, label = "Same r = -0.02 (p-value = 0.91)", size=3, color = "#0072B2", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = -2.5, y = 1.9, label = "Opposite r = -0.03 (p-value = 0.87)", size=3, color = "#D55E00", hjust=0) + ###This line needs to be manually changed
  annotate("text", x = -2.5, y = 1.8, label = "Mixed r = -0.03 (p-value = 0.78)", size=3, color = "#CC79A7", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

### Make scatter plot individual
### add pearson correlation manually
pdf("BHIVE_expression_2insertions_1gene_bydirection_same.pdf")
ggplot(same_gn, aes(x=V2, y=V3, color=group)) + 
  scale_color_manual(values=c("#0072B2")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression of 2 Insertions into 1 Gene") + 
  xlab("Insertion 1 (closest to TSS) Expression") + 
  ylab("Insertion 2 (closest to TES) Expression") + 
  scale_x_continuous(limits=c(-2.5, 2.5)) +
  scale_y_continuous(limits=c(-2.5, 2)) +
  annotate("text", x = -2.5, y = 2, label = "Same r = -0.02 (p-value = 0.91)", size=3, color = "#0072B2", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

### Make scatter plot individual
### add pearson correlation manually
pdf("BHIVE_expression_2insertions_1gene_bydirection_opposite.pdf")
ggplot(opp_gn, aes(x=V2, y=V3, color=group)) + 
  scale_color_manual(values=c("#D55E00")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression of 2 Insertions into 1 Gene") + 
  xlab("Insertion 1 (closest to TSS) Expression") + 
  ylab("Insertion 2 (closest to TES) Expression") + 
  scale_x_continuous(limits=c(-2.5, 2.5)) +
  scale_y_continuous(limits=c(-2.5, 2)) +
  annotate("text", x = -2.5, y = 2, label = "Opposite r = -0.03 (p-value = 0.87)", size=3, color = "#D55E00", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()

### Make scatter plot individual
### add pearson correlation manually
pdf("BHIVE_expression_2insertions_1gene_bydirection_mixed.pdf")
ggplot(mix_gn, aes(x=V2, y=V3, color=group)) + 
  scale_color_manual(values=c("#CC79A7")) +
  geom_point(alpha=1/4) + 
  geom_smooth(method=lm, se=FALSE) +
  ggtitle("HIV Expression of 2 Insertions into 1 Gene") + 
  xlab("Insertion 1 (closest to TSS) Expression") + 
  ylab("Insertion 2 (closest to TES) Expression") + 
  scale_x_continuous(limits=c(-2.5, 2.5)) +
  scale_y_continuous(limits=c(-2.5, 2)) +
  annotate("text", x = -2.5, y = 2, label = "Mixed r = -0.03 (p-value = 0.78)", size=3, color = "#CC79A7", hjust=0) + ###This line needs to be manually changed
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()
