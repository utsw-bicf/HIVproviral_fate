#!/bin/bash

#SBATCH --job-name=CHIP_SE
#SBATCH --partition=super
#SBATCH --output=chip.%j.out
#SBATCH --error=chip.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load nextflow/0.31.0
module add  python/3.6.1-2-anaconda

nextflow run workflow/main.nf \
--reads '/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/raw_data/*fastq.gz' \
--designFile '/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/design_chipseq_SE.txt' \
--genome 'GRCh38' \
--pairedEnd 'false' \
--outDir 'workflow/output_011519' \
-resume
