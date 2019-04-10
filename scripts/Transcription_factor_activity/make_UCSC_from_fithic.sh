#!/bin/bash

#SBATCH --job-name=fithic2UCSC
#SBATCH --partition=super
#SBATCH --output=fithic2UCSC.%j.out
#SBATCH --error=fithic2UCSC.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### This is A bash script to convert Fit-Hi-C output into visualization input for UCSC's Genome Browser in 'interact' format.
path='/project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/fitHiC'
files='rep1_1000000  rep1_20000  rep1_500000   rep2_150000  rep2_40000 rep1_150000   rep1_40000  rep2_1000000  rep2_20000   rep2_500000'
for file in ${files}; do
/home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Transcription_factor_activity/visualize-UCSC.sh \
  ${path}/${file}/FitHiC.spline_pass1.significances.txt.gz \
  ${path}/UCSC_browser_files/${file}.txt \
  0.05
done
