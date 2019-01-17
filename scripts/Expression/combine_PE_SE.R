### This R script combines the PE and SE RNAseq data
### combines fpkm, logCPM, countTable

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/combine_PE_SE')
SEdir <- '/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_SE/'
PEdir <- '/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_PE/'

### combine fpkm
SEfpkm <- read.table(file=paste(SEdir, "countTable.fpkm.txt", sep=""), header=T, sep='\t')
PEfpkm <- read.table(file=paste(PEdir, "countTable.fpkm.txt", sep=""), header=T, sep='\t')

SPfpkm <- merge(SEfpkm, PEfpkm, by=c("ENSEMBL", "SYMBOL", "TYPE"), all=T)
SPfpkm[is.na(SPfpkm)] <- 0
write.table(SPfpkm,file=("countTable.fpkm.txt"), quote=F ,row.names=F ,sep='\t')

### combine logCPM
SElogCPM <- read.table(file=paste(SEdir, "countTable.logCPM.txt", sep=""), header=T, sep='\t')
PElogCPM <- read.table(file=paste(PEdir, "countTable.logCPM.txt", sep=""), header=T, sep='\t')

SPlogCPM <- merge(SElogCPM, PElogCPM, by=c("ENSEMBL", "SYMBOL", "TYPE"), all=T)
SPlogCPM[is.na(SPlogCPM)] <- 0
write.table(SPlogCPM,file=("countTable.logCPM.txt"), quote=F ,row.names=F ,sep='\t')

### combine countTable
SEct <- read.table(file=paste(SEdir, "countTable.txt", sep=""), header=T, sep='\t')
PEct <- read.table(file=paste(PEdir, "countTable.txt", sep=""), header=T, sep='\t')

SPct <- merge(SEct, PEct, by=c("ENSEMBL", "SYMBOL", "TYPE"), all=T)
SPct[is.na(SPct)] <- 0
write.table(SPct,file=("countTable.txt"), quote=F ,row.names=F ,sep='\t')

### combine alignment summary
SEas <- read.table(file=paste(SEdir, "alignment.summary.txt", sep=""), header=T, sep='\t')
PEas <- read.table(file=paste(PEdir, "alignment.summary.txt", sep=""), header=T, sep='\t')

SPas <- rbind(SEas, PEas)
write.table(SPas,file=("alignment.summary.txt"), quote=F ,row.names=F ,sep='\t')
