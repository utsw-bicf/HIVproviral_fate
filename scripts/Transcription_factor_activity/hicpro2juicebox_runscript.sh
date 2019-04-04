#!/bin/bash

#SBATCH --job-name=h2j
#SBATCH --partition=super
#SBATCH --output=h2j.%j.out
#SBATCH --error=h2j.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load juicebox/1.5.6

/home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Transcription_factor_activity/hicpro2juicebox.sh \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/output_031419/hic_results/data/rep1/rep1.allValidPairs \
  -g /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/chrom.sizes \
  -j /cm/shared/apps/juicebox/1.5.6/Juicebox.jar



