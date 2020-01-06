#!/bin/bash

#SBATCH --job-name=cipher
#SBATCH --partition=super
#SBATCH --output=cipher.%j.out
#SBATCH --error=cipher.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load bedops/2.4.14


### This script uses the last part of the cipher enhancer identification package
### Rewritten in bash, just for this project

### Need DNase, H3K27Ac, H3K4me1, and H3K4me3

### Convert significant (-log10(0.05)=1.30102999566), replicated narrow peaks into beds
DNase='/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/consensusPeaks/DNAse-seq_ENCODE_GSM736501.replicated.narrowPeak'
H3K27Ac='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/consensusPeaks/H3K27ac.replicated.narrowPeak'
H3K4me1='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/GSM1603225_GFP_H3K4me1.replicated.narrowPeak'
H3K4me3='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/consensusPeaks/H3K4me3.replicated.narrowPeak'

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/cipher/bedfiles
awk '{OFS="\t"}; $9 > 1.30102999566 && $1 !~ /_/ {print $1,$2,$3,".",$4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/consensusPeaks/DNAse-seq_ENCODE_GSM736501.replicated.narrowPeak | sort-bed - > bedfiles/DNase.bed

awk '{OFS="\t"}; $9 > 1.30102999566 && $1 !~ /_/ {print $1,$2,$3,".",$4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/consensusPeaks/H3K27ac.replicated.narrowPeak | sort-bed - > bedfiles/H3K27Ac.bed

awk '{OFS="\t"}; $9 > 1.30102999566 && $1 !~ /_/ {print $1,$2,$3,".",$4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/GSM1603225_GFP_H3K4me1.replicated.narrowPeak | sort-bed - > bedfiles/H3K4me1.bed

awk '{OFS="\t"}; $9 > 1.30102999566 && $1 !~ /_/ {print $1,$2,$3,".",$4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/consensusPeaks/H3K4me3.replicated.narrowPeak | sort-bed - > bedfiles/H3K4me3.bed

###
