#!/bin/bash

#SBATCH --job-name=h2f
#SBATCH --partition=super
#SBATCH --output=h2f.%j.out
#SBATCH --error=h2f.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load python/2.7.5
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiCpro2fitHiC/rep1
cd /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiCpro2fitHiC/rep1
python /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Transcription_factor_activity/hicpro2fithic.py \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/output_031419/hic_results/matrix/rep1/raw/1000000/rep1_1000000.matrix \
  -b /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/output_031419/hic_results/matrix/rep1/raw/1000000/rep1_1000000_abs.bed \
  -s /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/output_031419/hic_results/matrix/rep1/iced/1000000/rep1_1000000_iced.matrix.biases

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiCpro2fitHiC/rep2
cd /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiCpro2fitHiC/rep2
python /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Transcription_factor_activity/hicpro2fithic.py \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/output_031419/hic_results/matrix/rep2/raw/1000000/rep2_1000000.matrix \
  -b /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/output_031419/hic_results/matrix/rep2/raw/1000000/rep2_1000000_abs.bed \
  -s /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/output_031419/hic_results/matrix/rep2/iced/1000000/rep2_1000000_iced.matrix.biases
