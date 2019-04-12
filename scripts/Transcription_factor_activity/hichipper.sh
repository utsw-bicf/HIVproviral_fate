#!/bin/bash

#SBATCH --job-name=hichipper
#SBATCH --partition=super
#SBATCH --output=hichipper.%j.out
#SBATCH --error=hichipper.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### Use the program "HiChipper" on the hiChiP data processed by HiC-pro
### Downloaded the program by
### module load python/2.7.x-anaconda
### pip install --user hichipper

###### Dependencies
module load python/2.7.x-anaconda
module load bedtools/2.26.0
module load samtools/gcc/1.8
module load R/3.5.1-gccmkl
source ~/.bash_profile 

###### Make yaml file
#echo "peaks:
#  - /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/SRR6010201_YY1_pooled_peaks.narrowPeak
resfrags:
#  - /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/GRCh38_MbolI.bed
#hicpro_output:
#  - /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/output_031419" >/project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiChipper/config.yaml

#hichipper --out YY1 /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiChipper/config.yaml

######### Call peaks
echo "peaks:
  - COMBINED,ALL
resfrags:
  - /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/GRCh38_MbolI.bed
hicpro_output:
  - /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/output_031419" >/project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiChipper/config_peaks_all.yaml

hichipper --out YY1_peaks_all /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiChipper/config_peaks_all.yaml


echo "peaks:
  - COMBINED,SELF
resfrags:
  - /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/GRCh38_MbolI.bed
hicpro_output:
  - /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/output_031419" >/project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiChipper/config_peaks_self.yaml

hichipper --out YY1_peaks_self /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiChipper/config_peaks_self.yaml
