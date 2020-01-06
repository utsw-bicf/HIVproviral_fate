### This R script plots the avg chromHMM score vs whether it's in a tad/loop

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib')

### Plot tads
tad_wHIV <- read.table("lib.tad.2D_filtered_scored.bed", header=F, sep='\t')
tad_woHIV <- read.table("noHIV_in_scored_tad.bed", header=F, sep='\t')

### Add identifider
tad_wHIV$label <- "HIV inserted"
tad_woHIV$label <- "No HIV inserted"

### Keep what maters and merge
tad_wHIV <- tad_wHIV[,c("label", "V5")]
tad_woHIV <- tad_woHIV[,c("label", "V5")]

df <- rbind(tad_wHIV, tad_woHIV)

### Plot
library(ggplot2)

pdf("HIV_insertion_avg_chromHMM_score_tad.pdf")
ggplot(df, aes(x=label, y=V5)) + 
  geom_violin(trim = FALSE) +
  ggtitle("HIV insertion, average chromHMM score in TADs") + 
  ylab("Average chromHMM score") +
  xlab("HIV status") +
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.5, color="red") +
  theme(plot.title = element_text(size=22, hjust = 0.5), 
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")
dev.off()