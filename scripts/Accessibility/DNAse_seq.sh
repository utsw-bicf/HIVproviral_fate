#!/bin/bash

#SBATCH --job-name=DN
#SBATCH --partition=super
#SBATCH --output=DN.%j.out
#SBATCH --error=DN.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load nextflow
module add  python/3.6.1-2-anaconda

### This pipeline is the atac-seq except
### --shift -100 --extsize 200 was hard coded into the call_peaks_macs.py
nextflow run workflow/main.nf \
  --reads '/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/raw_data/*.gz' \
  --designFile '/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/Design_DNAse_SE.txt' \
  --genome 'GRCh38' \
  --pairedEnd false \
  --tn5 false 
