### This R script plots a histogram of
### 1) length of ttseq vs FandR transcribed
### 2) Counts of ttseq vs FandR transcribed in 12,500 bp region

library(ggplot2)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/')

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
