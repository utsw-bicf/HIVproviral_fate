### This script makes QC plots for Hi-C data

library(ggplot2)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1/QCplots')

####################################################################################################
### Plot petagLocalDistribution
# Lib1
pLD <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1/petag.LocalDistribution.txt", header=T, sep="\t")

pdf("PE_tag_distance_lib1.pdf")
ggplot(pLD, aes(Local.Distance.in.bp.between.PE.tags)) + 
  geom_line(aes(y = Same.Strand, color="Same Strand")) + 
  geom_line(aes(y = Opposite.Strands, color="Opposite Strand")) +
  ylab("Read Counts") +
  xlab("Distance between HiC Read Pairs - lib1") +
  scale_x_continuous(limits = c(-1500, 1500)) +
  scale_color_discrete(name = "") #This removes legen title
dev.off()  

# Lib2
pLD <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2/petag.LocalDistribution.txt", header=T, sep="\t")

pdf("PE_tag_distance_lib2.pdf")
ggplot(pLD, aes(Local.Distance.in.bp.between.PE.tags)) + 
  geom_line(aes(y = Same.Strand, color="Same Strand")) + 
  geom_line(aes(y = Opposite.Strands, color="Opposite Strand")) +
  ylab("Read Counts") +
  xlab("Distance between HiC Read Pairs - lib2") +
  scale_x_continuous(limits = c(-1500, 1500)) +
  scale_color_discrete(name = "") #This removes legen title
dev.off()  

####################################################################################################
### Plot Restriction site distribution
# Lib1
pRS <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1/petagRestrictionDistribution.GATC.mis0.txt", header=T, sep="\t")

pdf("Restriction_site_distribution_lib1.pdf")
ggplot(pRS, aes(Distance.from.GATC.site)) + 
  geom_line(aes(y = Reads...strand, color="Reads + Strand")) + 
  geom_line(aes(y = Reads...strand.1, color="Reads - Strand")) +
  ylab("Read Counts - lib1") +
  xlab("Distance from GATC site") +
  scale_x_continuous(limits = c(-800, 800)) +
  scale_color_discrete(name = "") #This removes legen title
dev.off()  

# Lib2
pRS <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2/petagRestrictionDistribution.GATC.mis0.txt", header=T, sep="\t")

pdf("Restriction_site_distribution_lib2.pdf")
ggplot(pRS, aes(Distance.from.GATC.site)) + 
  geom_line(aes(y = Reads...strand, color="Reads + Strand")) + 
  geom_line(aes(y = Reads...strand.1, color="Reads - Strand")) +
  ylab("Read Counts - lib2") +
  xlab("Distance from GATC site") +
  scale_x_continuous(limits = c(-800, 800)) +
  scale_color_discrete(name = "") #This removes legen title
dev.off() 



####################################################################################################
### Plot petagFreqDistribution
pFD <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1/petag.FreqDistribution_1000.txt", header=T, sep="\t")

colnames(pFD) <- c("Distance", "Frequency")
pFD$Distance <- as.numeric(as.character(pFD$Distance))

pdf("PE_tag_distDist_lib1.pdf")
ggplot(data=pFD, aes(x=Distance, y=Frequency)) +
  geom_line() +
  scale_x_continuous(limits = c(0, 100000))
dev.off()

pFD <- read.table(file="/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2/petag.FreqDistribution_1000.txt", header=T, sep="\t")

colnames(pFD) <- c("Distance", "Frequency")
pFD$Distance <- as.numeric(as.character(pFD$Distance))

pdf("PE_tag_distDist_lib2.pdf")
ggplot(data=pFD, aes(x=Distance, y=Frequency)) +
  geom_line() +
  scale_x_continuous(limits = c(0, 100000))
dev.off()