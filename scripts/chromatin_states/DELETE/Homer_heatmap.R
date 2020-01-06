### From http://crazyhottommy.blogspot.com/2013/08/how-to-make-heatmap-based-on-chip-seq.html
library(gplots)
library(pheatmap)

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/homer_heatmaps')

d1 <- read.table("DNase_output_homer_TE.txt", header=T)

m1<- as.matrix( d1[,2:ncol(d1)])  # the first column is the TSS id,
rownames(m1)<- d1$Gene  # heatmap.2 works only on matrix, turn the dataframe to matrix, and                                                        # add the TSS id as the row name
m1<- log2(m1+1)  # log2 transform the raw counts

m.row.sum<- cbind(m1, rowSums(m1))
dim(m.row.sum) ### [1] 767  82
o1<- rev(order(m.row.sum[,82]))

m.row.sum<- m.row.sum[o1,]

bk = unique(c(seq(-1,3, length=100),seq(3,160,length=100)))

hmcols<- colorRampPalette(c("white","red"))(length(bk)-1)

pheatmap( m.row.sum[,1:82], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)
