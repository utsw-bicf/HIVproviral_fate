### Try calling enhancers on TTseq data using groHMM
### Use R 3.5, BiocManager::install("groHMM")
library(groHMM)
setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/groHMM')

tt1 <- as(readGAlignments(system.file("TTseq1", 
                                      "/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam", 
                                      package="groHMM")), "GRanges")