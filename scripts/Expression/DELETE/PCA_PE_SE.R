### Runs part of the dea analysis for RNAseq
### makes PCA plot

library(DESeq2)
library("RColorBrewer")
library("gplots")

####################
rowMax <- function(x) apply(x,1,max)
col.grp <- function(n,b) {
  colorvec <- vector(mode="character", length=length(n))
  plotcols <- rainbow(length(b))
  for (i in 1:length(n)) {
    for (j in 1:length(b)) {
      if ( n[i] == b[j] ) {
        colorvec[i] = plotcols[j]
      }
    }
  }
  c(colorvec)
}

####################

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/combine_PE_SE')

PEd <- read.table('/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/design_PE.txt', header=T, sep='\t')
SEd <- read.table('/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/design_SE.txt', header=T, sep='\t')
dtbl <- rbind(SEd, PEd)

tbl <- read.table('countTable.txt',header=TRUE,sep="\t")
tbl2 <- read.table('countTable.logCPM.txt',header=TRUE,sep="\t")
ct <- tbl[,4:length(tbl)]
row.names(ct) <- tbl$ENSEMBL

samples<- names(ct)

#dtbl <- read.table('design.txt',header=TRUE,sep="\t")
samtbl <- merge(as.data.frame(samples),dtbl,by.x="samples",by.y="SampleID",all.x=TRUE,sort=FALSE)
grpnames <- levels(factor(as.character(samtbl$SampleGroup)))
samtbl$SampleGroup <- factor(samtbl$SampleGroup, levels=grpnames)

colData <- samtbl[c('SampleGroup','SubjectID')]
row.names(colData) <- samtbl$samples

dds <- DESeqDataSetFromMatrix(countData=ct,colData= colData,design= ~ SampleGroup)
dds <- dds[ rowMax(counts(dds)) > 30, ]
dds <- dds[ colSums(counts(dds)) > 1000000]

if (nrow(counts(dds)) < 1) {
  print(paste("Samples are filtered if there is < 1M reads.  There are less than no remaining sample(s) after this filter",sep=' '))
  q()
}

countTable <- counts(dds)
grps <- as.character(colData(dds)$SampleGroup)
col.blocks <-col.grp(grps,levels(factor(grps)))
libSizes <- as.vector(colSums(countTable))

logcpm <- tbl2[,4:length(tbl2)]
row.names(logcpm) <- tbl2$SYMBOL

tmp.tab <- aggregate(t(logcpm),by=list(as.character(samtbl$SampleGroup)),FUN=mean)
row.names(tmp.tab) <- tmp.tab$Group.1
mean.by.group <- round(log2(t(tmp.tab[,c(2:ncol(tmp.tab))])+1), digits = 2)

#################### Run DESEQ2 ################################
dds <- DESeq(dds)
rld <- rlogTransformation(dds, blind=TRUE)
sampleDists <- dist(t(assay(rld)))

png(file="samples_heatmap.png",bg ="white",height=768,width=1024)
par(mar=c(7,4,4,2)+0.1) 
heatmap.2(as.matrix(sampleDists), col = bluered(100),RowSideColors = col.blocks,srtRow=45,srtCol=45,trace="none", margins=c(8,8), cexRow = 1.5, cexCol = 1.5)
dev.off()

#Compare Samples using PCA
png(file="pca.png",bg ="white",height=768,width=1024)
print(plotPCA(rld, intgroup="SampleGroup"),col.hab=col.blocks)
dev.off()