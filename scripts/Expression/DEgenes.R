### This script identifies genes that are DE between
### Emily and SRR213027*

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/combine_PE_SE')

eR <- read.table('Emily_SRR213027.edgeR.txt',header=TRUE,sep="\t")

### Filter for an fdr <- 0.05 and abs(logFC) >=1
filter1 <- eR[eR$fdr <= 0.05,]
filter2 <- filter1[abs(filter1$logFC) >= 1,]

### 8443 out of 15667 are Differentially expressed 