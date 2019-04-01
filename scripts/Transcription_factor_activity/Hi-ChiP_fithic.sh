#!/bin/bash

#SBATCH --job-name=fithic
#SBATCH --partition=super
#SBATCH --output=fithic.%j.out
#SBATCH --error=fithic.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load python/3.6.4-anaconda

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/fitHiC/rep1
python /home2/s185797/.local/lib/python3.6/site-packages/fithic/fithic.py \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiCpro2fitHiC/rep1/fithic.interactionCounts.gz \
  -f /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiCpro2fitHiC/rep1/fithic.fragmentMappability.gz \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/fitHiC/rep1 \
  -r 0 \
  -t /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiCpro2fitHiC/rep1/fithic.biases.gz \

  
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/fitHiC/rep2
python /home2/s185797/.local/lib/python3.6/site-packages/fithic/fithic.py \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiCpro2fitHiC/rep2/fithic.interactionCounts.gz \
  -f /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiCpro2fitHiC/rep2/fithic.fragmentMappability.gz \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/fitHiC/rep2 \
  -r 0 \
  -t /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiCpro2fitHiC/rep2/fithic.biases.gz \







##### General installation notes
#### module load python/3.6.4-anaconda
#### pip install --user fithic



