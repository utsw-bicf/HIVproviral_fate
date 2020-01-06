#!/bin/bash

#SBATCH --job-name=groHMM
#SBATCH --partition=super
#SBATCH --output=gHMM.%j.out
#SBATCH --error=gHMM.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load R/3.5.1-gccmkl

Rscript /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/chromatin_states/groHMM_venkat_wrapper.R --replicate1 /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam --replicate2 /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.dedup.bam --genome hg38 --out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/groHMM --ltprob -200 --uts 5
