########## This R script plots the expression of HIV in and out of Tads and Loops
########## And calculates probabilities that this isn't random

library(ggplot2)
library(scales)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_TadLoop')

### Load files
HT <- read.table("HIVexpression_inTad.bed", header=F, sep="\t")
HNT <- read.table("HIVexpression_NOTinTad.bed", header=F, sep="\t")

HL <- read.table("HIVexpression_inLoop.bed", header=F, sep="\t")
HNL <- read.table("HIVexpression_NOTinLoop.bed", header=F, sep="\t")

### Add Tad and Loop column
HT$compartment <- "Tad"
HNT$compartment <- "Not in Tad"

HL$compartment <- "Loop"
HNL$compartment <- "Not in Loop"

### Merge and plot
DF <- rbind(HT,HNT,HL,HNL)


DF$compartment <- factor(DF$compartment, levels = c("Tad","Not in Tad","Loop","Not in Loop"))
pdf("BHIVE_exp_by_TadLoop.pdf")
ggplot(DF, aes(x=compartment, y=V5)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(position=position_jitter(width=.2, height=0), size=0.5) +
  ggtitle("HIV Expression by Compartment") + 
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
#### Hypergeometric test - Loops
library(IRanges)
library(broom)

# Total counts
loop <- read.table("Loop.bed", header=F, sep="\t")
loop <- loop[,c("V1","V2","V3")]
loop$compartment <- "Loop"
Nloop <- read.table("not_loop.bed", header=F, sep="\t")
Nloop$compartment <- "Not in Loop"
Aloop <- rbind(loop,Nloop)

# Counts that we need to measure
Hloops <- rbind(HL, HNL)


total <- sum(Aloop$V3-Aloop$V2+1) # Calculate total
pr2 <- data.frame() # empty dataframe
for (i in c("Loop","Not in Loop")) {
  Hloopsr <- Hloops[Hloops$compartment == i,] # get compartment
  Aloopr <- Aloop[which(Aloop$compartment == i),]
  Alooptotal <- sum(Aloopr$V3-Aloopr$V2+1)
  NROW(Hloopsr)
  NROW(Hloops)
  Alooptotal/total
  hgt <- as.data.frame(t(rbind(phyper(NROW(Hloopsr), Alooptotal, (total-Alooptotal), NROW(Hloops), lower.tail=F), dhyper(NROW(Hloopsr), Alooptotal, (total-Alooptotal), NROW(Hloops), log=T))))
  hgt$state <- i
  hgt$actual <- NROW(Hloopsr)
  hgt$percentgenome <- format(Alooptotal/total,5)
  hgt$perHMM <- (NROW(Hloopsr))/(NROW(Hloops))
  pr2 <- rbind(pr2,hgt)
}

colnames(pr2) <- c("phyper", "log(dhyper)", "state", "insertion_count", "state_pct_genome", "pct_insertions")
write.table(pr2, file=("Loop_hypergeometric_test.txt"), quote=F, row.names = F, sep='\t')

######################################################################

########## Loop through states, chance
total <- sum(Aloop$V3-Aloop$V2+1) # Calculate total
pr1 <- data.frame()
for (i in c("Loop","Not in Loop")) {
  Hloopsr <- Hloops[Hloops$compartment == i,] # get compartment
  Aloopr <- Aloop[which(Aloop$compartment == i),]
  Alooptotal <- sum(Aloopr$V3-Aloopr$V2+1)
  NROW(Hloopsr)
  NROW(Hloops)
  Alooptotal/total
  pt <- tidy(prop.test(x=c((NROW(Hloopsr)),Alooptotal),n=c((NROW(Hloops)),total),p=NULL,alternative="two.sided", conf.level=0.95, correct=F))
  pt$state <- i
  pt$actual <- NROW(Hloopsr)
  pt$percentgenome <- format(Alooptotal/total,5)
  pt$perloop <- (NROW(Hloopsr))/(NROW(Hloops))
  pr1 <- rbind(pr1,pt)
}
rownames(pr1) <-
  c("Loop","Not in Loop")

##########
########## Loop through states, more
total <- sum(Aloop$V3-Aloop$V2+1) # Calculate total
prm <- data.frame()
for (i in c("Loop","Not in Loop")) {
  Hloopsr <- Hloops[Hloops$compartment == i,] # get compartment
  Aloopr <- Aloop[which(Aloop$compartment == i),]
  Alooptotal <- sum(Aloopr$V3-Aloopr$V2+1)
  NROW(Hloopsr)
  NROW(Hloops)
  Alooptotal/total
  pt <- tidy(prop.test(x=c((NROW(Hloopsr)),Alooptotal),n=c((NROW(Hloops)),total),p=NULL,alternative="greater", conf.level=0.95, correct=F))
  pt$state <- i
  pt$actual <- NROW(Hloopsr)
  pt$percentgenome <- format(Alooptotal/total,5)
  pt$perloop <- (NROW(Hloopsr))/(NROW(Hloops))
  prm <- rbind(prm,pt)
}
rownames(prm) <-
  c("Loop","Not in Loop")


##########
########## Loop through states, less
total <- sum(Aloop$V3-Aloop$V2+1) # Calculate total
prl <- data.frame()
for (i in c("Loop","Not in Loop")) {
  Hloopsr <- Hloops[Hloops$compartment == i,] # get compartment
  Aloopr <- Aloop[which(Aloop$compartment == i),]
  Alooptotal <- sum(Aloopr$V3-Aloopr$V2+1)
  NROW(Hloopsr)
  NROW(Hloops)
  Alooptotal/total
  pt <- tidy(prop.test(x=c((NROW(Hloopsr)),Alooptotal),n=c((NROW(Hloops)),total),p=NULL,alternative="less", conf.level=0.95, correct=F))
  pt$state <- i
  pt$actual <- NROW(Hloopsr)
  pt$percentgenome <- format(Alooptotal/total,5)
  pt$perloop <- (NROW(Hloopsr))/(NROW(Hloops))
  prl <- rbind(prl,pt)
}
rownames(prl) <-
  c("Loop","Not in Loop")

write.table(pr1, file=("Loop_prob.test_two.ways.txt"), quote=F, row.names = F, sep='\t')
write.table(prl, file=("Loop_prob.test_less.txt"), quote=F, row.names = F, sep='\t')
write.table(prm, file=("Loop_prob.test_more.txt"), quote=F, row.names = F, sep='\t')


### Check if mean expression is different; 2 way annova
# combine and reduce table
DF <- rbind(HL, HNL)
DFr <- DF[c("V5", "compartment")]

DFr$compartment <- factor(DFr$compartment,levels = c("Loop","Not in Loop"))
levels(DFr$compartment)

# Calculate summary stats
library(dplyr)
sumStats <- group_by(DFr, compartment) %>%
  summarise(
    count = n(),
    mean = mean(V5, na.rm = TRUE),
    sd = sd(V5, na.rm = TRUE)
  )

# Compute the analysis of variance
res.aov <- aov(V5 ~ compartment, data = DFr)

# Summary of the analysis
summary(res.aov)

### This result means that there ARE differences in the group
### Need to do a two-way to determine which group

# Tukey multiple pairwise-comparisons
TukeyHSD(res.aov)

# Multiple comparisons using multcomp package
library(multcomp)
summary(glht(res.aov, linfct = mcp(compartment = "Tukey")))

# Pairewise t-test
pairwise.t.test(DFr$V5, DFr$compartment, p.adjust.method = "BH")

# Check the homogeneity of variance assumption
library(car)
leveneTest(V5 ~ compartment, data = DFr)

# Check the normality assumption
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals )

### Data meets variance assumption but not normality
kruskal.test(V5 ~ compartment, data = DFr)

test <- pairwise.wilcox.test(DFr$V5, DFr$compartment, p.adjust.method = "BH")

test2 <- as.table(test$p.value)
write.table(test2, file="Loop_pairwise_wilcox_test.txt", quote=F, col.names = T, row.names = T, sep = "\t")









######################################################################
######################################################################
######################################################################
######################################################################
######################################################################

########## Tad
# Total counts
tad <- read.table("Tad.bed", header=F, sep="\t")
tad <- tad[,c("V1","V2","V3")]
tad$compartment <- "Tad"
Ntad <- read.table("not_tad.bed", header=F, sep="\t")
Ntad$compartment <- "Not in Tad"
Atad <- rbind(tad,Ntad)

# Counts that we need to measure
Htads <- rbind(HT, HNT)


total <- sum(Atad$V3-Atad$V2+1) # Calculate total
pr2 <- data.frame() # empty dataframe
for (i in c("Tad","Not in Tad")) {
  Htadsr <- Htads[Htads$compartment == i,] # get compartment
  Atadr <- Atad[which(Atad$compartment == i),]
  Atadtotal <- sum(Atadr$V3-Atadr$V2+1)
  NROW(Htadsr)
  NROW(Htads)
  Atadtotal/total
  hgt <- as.data.frame(t(rbind(phyper(NROW(Htadsr), Atadtotal, (total-Atadtotal), NROW(Htads), lower.tail=F), dhyper(NROW(Htadsr), Atadtotal, (total-Atadtotal), NROW(Htads), log=T))))
  hgt$state <- i
  hgt$actual <- NROW(Htadsr)
  hgt$percentgenome <- format(Atadtotal/total,5)
  hgt$perHMM <- (NROW(Htadsr))/(NROW(Htads))
  pr2 <- rbind(pr2,hgt)
}

colnames(pr2) <- c("phyper", "log(dhyper)", "state", "insertion_count", "state_pct_genome", "pct_insertions")
write.table(pr2, file=("Tad_hypergeometric_test.txt"), quote=F, row.names = F, sep='\t')

######################################################################

########## Loop through states, chance
total <- sum(Atad$V3-Atad$V2+1) # Calculate total
pr1 <- data.frame()
for (i in c("Tad","Not in Tad")) {
  Htadsr <- Htads[Htads$compartment == i,] # get compartment
  Atadr <- Atad[which(Atad$compartment == i),]
  Atadtotal <- sum(Atadr$V3-Atadr$V2+1)
  NROW(Htadsr)
  NROW(Htads)
  Atadtotal/total
  pt <- tidy(prop.test(x=c((NROW(Htadsr)),Atadtotal),n=c((NROW(Htads)),total),p=NULL,alternative="two.sided", conf.level=0.95, correct=F))
  pt$state <- i
  pt$actual <- NROW(Htadsr)
  pt$percentgenome <- format(Atadtotal/total,5)
  pt$pertad <- (NROW(Htadsr))/(NROW(Htads))
  pr1 <- rbind(pr1,pt)
}
rownames(pr1) <-
  c("Tad","Not in Tad")

##########
########## Loop through states, more
total <- sum(Atad$V3-Atad$V2+1) # Calculate total
prm <- data.frame()
for (i in c("Tad","Not in Tad")) {
  Htadsr <- Htads[Htads$compartment == i,] # get compartment
  Atadr <- Atad[which(Atad$compartment == i),]
  Atadtotal <- sum(Atadr$V3-Atadr$V2+1)
  NROW(Htadsr)
  NROW(Htads)
  Atadtotal/total
  pt <- tidy(prop.test(x=c((NROW(Htadsr)),Atadtotal),n=c((NROW(Htads)),total),p=NULL,alternative="greater", conf.level=0.95, correct=F))
  pt$state <- i
  pt$actual <- NROW(Htadsr)
  pt$percentgenome <- format(Atadtotal/total,5)
  pt$pertad <- (NROW(Htadsr))/(NROW(Htads))
  prm <- rbind(prm,pt)
}
rownames(prm) <-
  c("Tad","Not in Tad")


##########
########## Loop through states, less
total <- sum(Atad$V3-Atad$V2+1) # Calculate total
prl <- data.frame()
for (i in c("Tad","Not in Tad")) {
  Htadsr <- Htads[Htads$compartment == i,] # get compartment
  Atadr <- Atad[which(Atad$compartment == i),]
  Atadtotal <- sum(Atadr$V3-Atadr$V2+1)
  NROW(Htadsr)
  NROW(Htads)
  Atadtotal/total
  pt <- tidy(prop.test(x=c((NROW(Htadsr)),Atadtotal),n=c((NROW(Htads)),total),p=NULL,alternative="less", conf.level=0.95, correct=F))
  pt$state <- i
  pt$actual <- NROW(Htadsr)
  pt$percentgenome <- format(Atadtotal/total,5)
  pt$pertad <- (NROW(Htadsr))/(NROW(Htads))
  prl <- rbind(prl,pt)
}
rownames(prl) <-
  c("Tad","Not in Tad")

write.table(pr1, file=("Tad_prob.test_two.ways.txt"), quote=F, row.names = F, sep='\t')
write.table(prl, file=("Tad_prob.test_less.txt"), quote=F, row.names = F, sep='\t')
write.table(prm, file=("Tad_prob.test_more.txt"), quote=F, row.names = F, sep='\t')


### Check if mean expression is different; 2 way annova
# combine and reduce table
DF <- rbind(HT, HNT)
DFr <- DF[c("V5", "compartment")]

DFr$compartment <- factor(DFr$compartment,levels = c("Tad","Not in Tad"))
levels(DFr$compartment)

# Calculate summary stats
library(dplyr)
sumStats <- group_by(DFr, compartment) %>%
  summarise(
    count = n(),
    mean = mean(V5, na.rm = TRUE),
    sd = sd(V5, na.rm = TRUE)
  )

# Compute the analysis of variance
res.aov <- aov(V5 ~ compartment, data = DFr)

# Summary of the analysis
summary(res.aov)

### This result means that there ARE differences in the group
### Need to do a two-way to determine which group

# Tukey multiple pairwise-comparisons
TukeyHSD(res.aov)

# Multiple comparisons using multcomp package
library(multcomp)
summary(glht(res.aov, linfct = mcp(compartment = "Tukey")))

# Pairewise t-test
pairwise.t.test(DFr$V5, DFr$compartment, p.adjust.method = "BH")

# Check the homogeneity of variance assumption
library(car)
leveneTest(V5 ~ compartment, data = DFr)

# Check the normality assumption
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals )

### Data meets variance assumption but not normality
kruskal.test(V5 ~ compartment, data = DFr)

test <- pairwise.wilcox.test(DFr$V5, DFr$compartment, p.adjust.method = "BH")

test2 <- as.table(test$p.value)
write.table(test2, file="Tad_pairwise_wilcox_test.txt", quote=F, col.names = T, row.names = T, sep = "\t")


