### This R script plots a histogram of
### 1) length of ttseq vs FandR transcribed
### 2) Counts of ttseq vs FandR transcribed in 12,500 bp region

library(ggplot2)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/')

enh <- read.table("enhancer_150min_merged12500_rpkm_filtered_peaks.tsv", header=T, sep="\t")

enh$fwd <- enh$ttfwd87+enh$ttfwd88
enh$rev <- enh$ttrev87+enh$ttrev88
enh$total <- enh$fwd+enh$rev
enh$length <- enh$end-enh$start
enh$diff <- round(abs(enh$fwd/enh$rev),digits=2)

hist(enh$total)
hist(enh$length)
enh2.b <- enh[which(enh$total > 1),]
hist(enh2.b$total, breaks=150)
enh2 <- enh[which(enh$length > 9000),]
hist(enh2$length)
hist(enh2$total)
ggplot(enh2, aes(x=enh2$length, y=enh2$total))

### Sort low to high on length and add rank
enh2S <- enh2[order(enh2$total),]
enh2S$rank <- 1:nrow(enh2S)

ggplot(enh2S, aes(x=rank, y=total)) +
  geom_point(shape=1) 


enh3 <- enh[which(enh$total > 1),]
enh4 <- enh3[which(enh3$length > 9000),]
enh4$bm <- round(enh4$total*enh4$diff,digits = 2)
enh4S <- enh4[order(enh4$total),]
enh4S$rank <- 1:nrow(enh4S)

ggplot(enh4S, aes(x=rank, y=total)) +
  geom_point(shape=1) 

# based on diff
enh4S <- enh4[order(enh4$diff),]
enh4S$rank <- 1:nrow(enh4S)

ggplot(enh4S, aes(x=rank, y=total)) +
  geom_point(shape=1) 

enh5 <- enh4S[which(enh4$total > 1 & enh4$diff <1),]
ggplot(enh4S, aes(x=total, y=diff)) +
  geom_point(shape=1) 
enh6 <- enh4S[which(enh4S$total > 12.5 & enh4S$diff < 3),]
ggplot(enh6, aes(x=total, y=diff)) +
  geom_point(shape=1) 




A <- read.table("enhancers_150bp.bed", header=F, sep="\t")

A$f <- A$V4+A$V5
A$r <- A$V6+A$V7
A$t <- A$f+A$r
A$l <- A$V3-A$V2+1


hist(A$t)
hist(A$l)
ggplot(A, aes(x=A$l, y=A$t)) + geom_point(shape=1)

As <- A[which(A$l >9000),]
ggplot(As, aes(x=As$l, y=As$t)) + geom_point(shape=1)

ttrpkm2.b <- ttrpkm[which(ttrpkm$total > 1),]
hist(ttrpkm2.b$total, breaks=150)
ttrpkm2 <- ttrpkm[which(ttrpkm$length > 9000),]
hist(ttrpkm2$length)
hist(ttrpkm2$total)
ggplot(ttrpkm2, aes(x=ttrpkm2$length, y=ttrpkm2$total))

### Sort low to high on length and add rank
ttrpkm2S <- ttrpkm2[order(ttrpkm2$length),]
ttrpkm2S$rank <- 1:nrow(ttrpkm2S)

ggplot(ttrpkm2S, aes(x=rank, y=total)) +
  geom_point(shape=1)


######################Not 12.5
ori <- read.table("jurkat_2_notgenes_transcribed_inFandR.bed", header=T, sep="\t")
ori$length <- ori$X34200-ori$X11600+1
######################Not 12.5


######################Annotations
ann <- read.table("jurkat_2_notgenes_transcribed_inFandR_annotated2.bed", header=F, sep="\t")
count <- table(ann$V15, ann$V20)
cols = c("dodgerblue2", "#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black", "gold1",
         "skyblue2", "#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
         "darkturquoise", "green1", "yellow4", "yellow3",
         "darkorange4", "brown")
pdf("test.pdf")
barplot(count, col= cols)
dev.off()
pdf("test2.pdf")
barplot(count, col= cols, legend.text = T)
dev.off()

### Make a stacked bargraph of chromHMM with type



######################Annotations


###
ttrpkm <- read.table("jurkat_2_notgenes_transcribed_inFandR_merged12500_rpkm_filtered_peaks.tsv", header=T, sep="\t")
ttrpkm$fwd <- ttrpkm$ttfwd87+ttrpkm$ttfwd88
ttrpkm$rev <- ttrpkm$ttrev87+ttrpkm$ttrev88
ttrpkm$total <- ttrpkm$fwd+ttrpkm$rev
ttrpkm$length <- ttrpkm$end-ttrpkm$start

hist(ttrpkm$total)
hist(ttrpkm$length)
ttrpkm2.b <- ttrpkm[which(ttrpkm$total > 1),]
hist(ttrpkm2.b$total, breaks=150)
ttrpkm2 <- ttrpkm[which(ttrpkm$length > 9000),]
hist(ttrpkm2$length)
hist(ttrpkm2$total)
ggplot(ttrpkm2, aes(x=ttrpkm2$length, y=ttrpkm2$total))

### Sort low to high on length and add rank
ttrpkm2S <- ttrpkm2[order(ttrpkm2$total),]
ttrpkm2S$rank <- 1:nrow(ttrpkm2S)

ggplot(ttrpkm2S, aes(x=rank, y=total)) +
  geom_point(shape=1) 




ttrpkm2S <- ttrpkm2[order(ttrpkm2$length),]
ttrpkm2S$rank <- 1:nrow(ttrpkm2S)

ggplot(ttrpkm2S, aes(x=rank, y=total)) +
  geom_point(shape=1)






ggplot(ttrpkm2, aes(x=length, y=total)) +
  geom_point(shape=1)









### Load ttseq vs FandR transcribed file
tt1 <- read.table("jurkat_2_notgenes_transcribed_inFandR.bed", header=F, sep="\t")
### remove alternative chromosomes
tt1 <- tt1[!grepl("_",tt1$V1),]
### remove identical rows
tt1 <- unique(tt1[ , 1:3 ])

tt1$diff <- tt1$V3-tt1$V2+1

hist(tt1$diff, breaks=20)


### Load count data
ttc <- read.table("jurkat_2_notgenes_transcribed_inFandR.countsin12500.bed", header=F, sep="\t")
### remove alternative chromosomes
ttc <- ttc[!grepl("_",ttc$V1),]
### remove counts of 0
ttc <- ttc[!(ttc$V4 == 0),]

### Subset for V4 >3
ttc2 <- ttc[(ttc$V4 > 6),]
hist(ttc2$V4)

### Remove counts that are overlapping then plot
library(data.table)
library(IRanges)
DT <- as.data.table(ttc)

## find interval for each chromsom
test <- DT[,group := { 
  ir <-  IRanges(V2, V3);
  subjectHits(findOverlaps(ir, reduce(ir)))
},by=V1]

test2 <- unique(test[,4:5])
hist(test2$V4)


### extend range by 12,500
tt1$start <- tt1$V2-12500
tt1$end <- tt1$V3+12500
# find overlapping regions
DT <- as.data.table(tt1)
itest <- DT[,group := { 
  ir <-  IRanges(start, end);
  subjectHits(findOverlaps(ir, reduce(ir)))
},by=V1]


is.merged.region(tt1)
btest <- bedr.merge.region(tt1,distance=12500)

itest2 <- itest[,(unique(itest[,c(1,7)]))]
