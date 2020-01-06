### This will run through the dea script
### Find DEgene differences between Emily and SRR213027*
### They are the same cell type and we shouldn't see much of a difference

#!/cm/shared/apps/R/intel/3.2.1/bin/Rscript

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/combine_PE_SE')

library(edgeR)
library(DESeq2)
library("RColorBrewer")
library("gplots")
library(qusage)

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


#################### Read in Data ################################
genenames <- read.table(file="/project/shared/bicf_workflow_ref/GRCh38/genenames.txt",header=TRUE,sep='\t')

tbl <- read.table('countTable.txt',header=TRUE,sep="\t")
tbl <- tbl[,c("ENSEMBL", "SYMBOL", "TYPE", "SRR2130272", "SRR2130273", "SRR2130274", "Emily_jurkat_A", "Emily_jurkat_B", "Emily_jurkat_C")]
tbl2 <- read.table('countTable.logCPM.txt',header=TRUE,sep="\t")
tbl2 <- tbl2[,c("ENSEMBL", "SYMBOL", "TYPE", "SRR2130272", "SRR2130273", "SRR2130274", "Emily_jurkat_A", "Emily_jurkat_B", "Emily_jurkat_C")]
ct <- tbl[,4:length(tbl)]
row.names(ct) <- tbl$ENSEMBL

samples<- names(ct)

### design table
PEd <- read.table('/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/design_PE.txt', header=T, sep='\t')
SEd <- read.table('/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/design_SE.txt', header=T, sep='\t')
ori_dtbl <- rbind(SEd, PEd)
write.table(ori_dtbl,file=("design_PE_SE.txt"), quote=F ,row.names=F ,sep='\t')
### Manually modify the output and use it further
dtbl <- read.table('design_PE_SE_4edgeR.txt',header=TRUE,sep="\t")

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

###################### Run EdgeR ########################
###################### Run QuSage ######################

design <- model.matrix(~grps)
d <- DGEList(counts=countTable,group=grps,lib.size=libSizes)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
png(file="mds.png",bg ="white",height=768,width=1024)
plotMDS(d, labels=grps,col=col.blocks, cex.axis=1.5, cex.lab=1.5, cex=1.5)
op <- par(cex = 1.5)
legend("topleft",legend=grpnames,col=rainbow(length(grpnames)),pch=20)
dev.off()
cond <-levels(d$samples$group)
colnames(design) <- levels(d$samples$group)
a <- length(cond)-1
for (i in 1:a) {
  for (j in 2:length(cond)) {
    if (i == j) {
      next
    } else {
      c <- exactTest(d, pair=c(cond[j],cond[i]))
      res <- c$table
      res2 <- merge(genenames,res,by.x='ensembl',by.y='row.names',all.y=TRUE,all.x=FALSE)
      output <- merge(res2,mean.by.group,by.y="row.names",by.x='symbol')
      output$rawP <- output$PValue
      output$logFC <- output$logFC
      output$fdr <- p.adjust(output$rawP, method ='fdr')
      output$bonf <- p.adjust(output$rawP, method ='bonferroni')
      write.table(output,file=paste(cond[i],'_',cond[j],'.edgeR.txt',sep=""),quote=FALSE,row.names=FALSE,sep='\t')
      filt.out <- na.omit(output[output$fdr < 0.05,])
      if (nrow(filt.out) > 2) {
        subset <- logcpm[row.names(logcpm) %in% filt.out$symbol,]
        subset <- subset[!apply(subset, 1, function(x) {any(x == 0)}),]
        gnames <- filt.out[c('ensembl','symbol')]
        s <- merge(gnames,subset,by.x="ensembl",by.y="row.names",all.x=FALSE,all.y=TRUE,sort=FALSE)
        STREE <- hclust(dist(t(subset)))
        zscores <- scale(t(subset))
        ngenes <- length(colnames(zscores))
        textscale <- (1/(ngenes/30))
        if (textscale > 1) {
          textscale <-1
        }
        if (textscale < 0.1) {
          textscale <- 0.1
        }
        png(file=paste(cond[i],'_',cond[j],'.heatmap.edgeR.png',sep=""),height=768,width=1024)
        heatmap.2(zscores, col = bluered(100),Rowv = as.dendrogram(STREE), RowSideColors = col.blocks,dendrogram='row', cexCol=textscale,labCol=s$symbol,srtRow=45,srtCol=45,trace="none", margins=c(5, 5))
        legend("topright",legend=grpnames,col=rainbow(length(grpnames)),pch=20,cex=0.5)
        dev.off()
      }
    }
  }
}




