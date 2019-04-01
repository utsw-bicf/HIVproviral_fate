#!/bin/bash

#SBATCH --job-name=chromHMM
#SBATCH --partition=super
#SBATCH --output=HMM2.%j.out
#SBATCH --error=HMM2.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

#### This script runs chromHMM on chipseq data

#### Unzip and move the bed files
#for testing: names='SRR1603650 SRR1603649 SRR1603654 SRR1603652 SRR577483 SRR577484 SRR061743 SRR061744'
#names='SRR2043614 SRR1603650 SRR1603654 SRR577482 SRR577483 GSM1603213 GSM1603225 GSM1603215 GSM1603211 GSM1603209 GSM1603227 GSM1603217 GSM1603221 GSM1603223 SRR2157609 SRR1522115 SRR3722566 SRR3722577 SRR6010201 SRR2043612 SRR1603649 SRR1603652 SRR577484 SRR577484 GSM1603229 SRR1057276 SRR1522118 SRR3722570 SRR3722573 SRR6010202'


#for name in ${names}; do
#zcat /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/convertReads/${name}.filt.nodup.bedse.gz >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bed/${name}.bed
#done

#### Run Chipseq only first; then add other data types
### First make binary file
#java -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar BinarizeBed \
#  /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bed \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/chipseq_input.txt \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output

#for i in $(seq 5 20); do
#java -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar LearnModel \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_learn${i} \
#  ${i} \
#  hg38 
#done


unset DISPLAY
for i in $(seq 5 20); do
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar LearnModel \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_learn${i} \
  ${i} \
  hg38 
done

### The above didn't seem to work correctly. Try with just the histone
###########################################################################################################################

