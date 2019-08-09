#!/bin/Rscript

# Load libraries
if(!require("argparse")){
    install.packages("argparse", repos="http://cran.r-project.org")
    library(argparse)
}
if(!require("jsonlite")){
    install.packages("jsonlite", repos="http://cran.r-project.org")
    library("jsonlite")
}
source("http://bioconductor.org/biocLite.R")
if (!require("groHMM")){
    biocLite("groHMM")
    library("groHMM")
}
if (!require("GenomicFeatures")){
    biocLite("GenomicFeatures")
    library("GenomicFeatures")
}
if (!require("rtracklayer")){
    biocLite("rtracklayer")
    library("rtracklayer")
}
#if(!require("jsonlite")){
#    install.packages("jsonlite", repos="http://cran.r-project.org")
#    library("jsonlite")
#}

# Currently mouse or human
if (!require("TxDb.Hsapiens.UCSC.hg38.knownGene")){
    biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")
    library("TxDb.Hsapiens.UCSC.hg38.knownGene")
}
if (!require("TxDb.Hsapiens.UCSC.hg19.knownGene")){
    biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
    library("TxDb.Hsapiens.UCSC.hg19.knownGene")
}
if (!require("TxDb.Mmusculus.UCSC.mm10.knownGene")){
    biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
    library("TxDb.Mmusculus.UCSC.mm10.knownGene")
}
if (!require("org.Hs.eg.db")){
    biocLite("org.Hs.eg.db")
    library("org.Hs.eg.db")
}
if (!require("org.Mm.eg.db")){
    biocLite("org.Mm.eg.db")
    library("org.Mm.eg.db")
}

# Create parser object
parser <- ArgumentParser()

# Specify our desired options
parser$add_argument("-r1", "--replicate1", nargs=1, help = "File paths to Replicate 1", required = TRUE)
parser$add_argument("-r2", "--replicate2" ,nargs=1, help = "File paths to Replicate 2", required = TRUE)
parser$add_argument("-g", "--genome", help = "The ucsc genome", required = TRUE)
parser$add_argument("-o", "--out", help = "The output directory", required = TRUE)
parser$add_argument("-b", "--ltprob", type="integer", help = "The LtsProbB", required = TRUE)
parser$add_argument("-u", "--uts", type="integer", help = "The UTS", required = TRUE)

# Parse arguments
args <- parser$parse_args()

# Set the working directory to output directory
setwd(args$out)

# Set mc.cores to 4
options(mc.cores=getCores(4))
options(mc.cores=getCores(2))


# Load alignment files
alignment_1 <- as(readGAlignments(args$replicate1), "GRanges")
alignment_2 <- as(readGAlignments(args$replicate2), "GRanges")
alignment_1 <- as(readGAlignments("/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam"), "GRanges")
alignment_2 <- as(readGAlignments("/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.dedup.bam"), "GRanges")
alignment_1 <- as(readGAlignments("/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam"), "GRanges")
alignment_2 <- as(readGAlignments("/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.dedup.bam"), "GRanges")
alignments <- sort(c(alignment_1, alignment_2))

# Load UCSC Genes
if(args$genome=='hg19') {
    kgdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    org <- org.Hs.eg.db
} else if(args$genome=='mm10')  {
    kgdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    org <- org.Mm.eg.db
} else if(args$genome=='hg38')  {
    kgdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    org <- org.Hs.eg.db
}
kgtx <- transcripts(kgdb,columns=c("gene_id", "tx_id", "tx_name"))

# Collapse overlapping annotations
kgConsensus <- makeConsensusAnnotations(kgtx, keytype="gene_id",
        mc.cores=getOption("mc.cores"))
map <- select(org, keys=unlist(mcols(kgConsensus)$gene_id),
        columns=c("SYMBOL"), keytype=c("ENTREZID"))
mcols(kgConsensus)$symbol <- map$SYMBOL
mcols(kgConsensus)$type <- "gene"

# Detect transcripts
hmmResult <- detectTranscripts(alignments, LtProbB=args$ltprob, UTS=args$uts)
hmmResult <- detectTranscripts(alignments)





#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
alignments <- sort(alignment_1)

testa <- alignments[seqnames(alignments) == "chr1" & start(alignments) > 10000 & end(alignments) < 1000000]

hmmResult <- detectTranscripts(testa)
Fp <- testa[strand(testa) == "+"]
Fm <- testa[strand(testa) == "-"]
hmmResult <- detectTranscripts(Fp=Fp, Fm=Fm, threshold = 1)


nFp <- NROW(Fp)
nFm <- NROW(Fm)
CHRp <- as.character(names(Fp))
CHRm <- as.character(names(Fm))
size <- as.integer(50)
ANS <- NULL
HMM <- list()
HMM$nstates <- as.integer(2)
HMM$ePrDist <- c("dgamma", "dgamma")
HMM$iProb <- as.double(log(c(1, 0)))
UTS=5
HMM$ePrVars <- as.list(data.frame(c(UTS, 1/UTS, -1), c(0.5, 
                                                       10, -1)))
LtProbA = -5
LtProbB = -200
HMM$tProb <- as.list(data.frame(c(log(1 - exp(LtProbA)), 
                                  LtProbA), c(LtProbB, log(1 - exp(LtProbB)))))
FT <- list()
#for (i in 1:nFp) FT[[i]] <- as.double(Fp[[i]] + 1)

for (i in 1:nFp) FT[i] <- as.double(Fp[i] + 1)

for (i in 1:nFm) FT[[i + nFp]] <- as.double(Fm[[i]] + 1)

remove(Fp)
remove(Fm)
BWem <- .Call("RBaumWelchEM", HMM$nstates, FT, as.integer(1), 
              HMM$ePrDist, HMM$ePrVars, HMM$tProb, HMM$iProb, as.double(threshold), 
              c(TRUE, FALSE), c(FALSE, TRUE), as.integer(1), TRUE, 
              PACKAGE = "groHMM")
for (i in seq_along(CHRp)) {
  ans <- .Call("getTranscriptPositions", as.double(BWem[[3]][[i]]), 
               as.double(0.5), size, PACKAGE = "groHMM")
  Nrep <- NROW(ans$Start)
  ANS <- rbind(ANS, data.frame(chrom = rep(CHRp[i], Nrep), 
                               start = ans$Start, end = ans$End, strand = rep("+", 
                                                                              Nrep)))
}
for (i in seq_along(CHRm)) {
  ans <- .Call("getTranscriptPositions", as.double(BWem[[3]][[i + 
                                                                nFp]]), as.double(0.5), size, PACKAGE = "groHMM")
  Nrep <- NROW(ans$Start)
  ANS <- rbind(ANS, data.frame(chrom = rep(CHRm[i], NROW(ans$Start)), 
                               start = ans$Start, end = ans$End, strand = rep("-", 
                                                                              Nrep)))
}
BWem[[4]] <- GRanges(seqnames = Rle(ANS$chrom), ranges = IRanges(ANS$start, 
                                                                 ANS$end - 1), strand = Rle(strand(ANS$strand)), type = Rle("tx", 
                                                                                                                            NROW(ANS)), ID = paste(ANS$chrom, "_", ANS$start, ANS$strand, 
                                                                                                                                                   sep = ""))
names(BWem) <- c("emisParams", "transParams", "viterbiStates", 
                 "transcripts")
if (debug) {
  print(BWem[[1]])
  print(BWem[[2]])
}
return(BWem)
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################











txHMM <- hmmResult$transcripts

getExpressedAnnotations <- function(features, reads) {
       fLimit <- limitToXkb(features)
       count <- countOverlaps(fLimit, reads)
       features <- features[count!=0,]
       return(features[(quantile(width(features), .05) < width(features))
      & (width(features) < quantile(width(features), .95)),])}

conExpressed <- getExpressedAnnotations(features=kgConsensus,reads=alignments)

# Plot Density
#png('transcript-density-plot.png') # TODO: Fix plot output
td <- getTxDensity(txHMM, conExpressed, mc.cores=getOption("mc.cores"), plot=FALSE)
#u <- par("usr")
#lines(c(u[1], 0, 0, 1000, 1000, u[2]), c(0,0,u[4]-.04,u[4]-.04,0,0),col="red")
#legend("topright", lty=c(1,1), col=c(2,1), c("ideal", "groHMM"))
#text(c(-500,500), c(.05,.5), c("FivePrimeFP", "TP"))
#dev.off()
td_JSON <- toJSON(td, null='null',auto_unbox = TRUE, pretty=TRUE)
write(td_JSON,"transcript-density-quality-metrics.json", ncolumns = 1, append = FALSE, sep = " ")

# Repairing called transcripts
break_plus <- breakTranscriptsOnGenes(txHMM, kgConsensus, strand="+")
break_minus <- breakTranscriptsOnGenes(txHMM, kgConsensus, strand="-")
txBroken <- c(break_plus, break_minus)
txFinal <- combineTranscripts(txBroken, kgConsensus)
export(txFinal, con="final-transcripts.bed")

