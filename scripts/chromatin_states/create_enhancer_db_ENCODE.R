### This script puts together the output from creat_enhancer_db_ENCODE_HOMER.sh
library(dplyr)
library(data.table)
library(ggplot2)

### Load the bed files
#dataFiles <- lapply(Sys.glob("/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/bed_counts/*.bed"), read.table)
file1 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/bed_counts/ENCFF031NSV.counts.bed", header=F, sep="\t")
file2 <- read.table("/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/bed_counts/ENCFF083VQO.counts.bed", header=F, sep="\t")
colnames(file1) <- c("chr","start","stop",mget(file1))
colnames(file2) <- c("chr","start","stop","file2")
out <-file1  %>%
  full_join(file2,by = c("chr","start","stop")) 

out <-file1  %>%
  full_join(file2,by = c("V1","V2","V3")) 

file.names <- dir(path, pattern =".txt")
for(i in 1:length(file.names)){
  file <- read.table(file.names[i],header=TRUE, sep=";", stringsAsFactors=FALSE)
  out.file <- rbind(out.file, file)
}

### Load files and merge into dataframe
setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/bed_counts')
out <- as.data.frame()
file.names <- dir("/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/bed_counts", pattern ="*.counts.bed")
for(i in 1:length(file.names)){
  file <- read.table(file.names[i],header=F, sep="\t", stringsAsFactors=F)
  out <- out %>%
    full_join(file,by = c("V1","V2","V3"))
}

### Remove columns that are identical (mistakes); 10=11 and 16=17
out2 <- out[-c(11,16)]
### Remove rows that have all zero counts
out3 <- filter(out2, V4 > 0 & V4.x >0 & V4.y >0 & 
                 V4.x.x >0 & V4.y.y & 
                 V4.x.x.x >0 & V4.y.y.y & 
                 V4.x.x.x.x >0 & 
                 V4.x.x.x.x.x >0 & V4.y.y.y.y.y &
                 V4.x.x.x.x.x.x >0 & V4.y.y.y.y.y.y &
                 V4.y.y.y.y.y.y.y &
                 V4.x.x.x.x.x.x.x.x >0 & V4.y.y.y.y.y.y.y.y &
                 V4.x.x.x.x.x.x.x.x.x >0 & V4.y.y.y.y.y.y.y.y.y &
                 V4.x.x.x.x.x.x.x.x.x.x >0 & V4.y.y.y.y.y.y.y.y.y.y &
                 V4.x.x.x.x.x.x.x.x.x.x.x >0 & V4.y.y.y.y.y.y.y.y.y.y.y &
                 V4.x.x.x.x.x.x.x.x.x.x.x.x >0 & V4.y.y.y.y.y.y.y.y.y.y.y.y)
  
### Make histogram to find cutoff
dfhist1 <- out3[-c(1:3)]
dfhist2 <- melt(dfhist1)
ggplot(dfhist2, aes(x = variable, y = value)) +
  geom_boxplot()

### Keep rows where 80% of the columns have a value > 4 (cutoff based on histogram)
### Filter logCPMtbl for at least > 1 in at least 2 replicates of 1 treatment
test <- apply(out3[4:26], 1, function(x) sum(x >=5))
keepdf <- out3[which(test >= 18),]

keepdfr <- keepdf[c(1:3)]
colnames(keepdfr) <- c("chr", "start", "end")
write.table(keepdfr, file=("unmerged_coordinates.bed"), quote=F, col.names = F, row.names = F, sep="\t")
