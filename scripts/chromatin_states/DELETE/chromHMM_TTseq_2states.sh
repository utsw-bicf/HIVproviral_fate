#!/bin/bash

#SBATCH --job-name=TTchromHMM
#SBATCH --partition=super
#SBATCH --output=ttHMM.%j.out
#SBATCH --error=ttHMM.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END


### TT-seq
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/Ttseq_CRAMER_GSM2260187.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/Ttseq_CRAMER_GSM2260188.bam


#mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/binary
#java -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar BinarizeBam \
#  /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/ttseq_2state_input_bam.txt \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/binary

### Run the modeling
#unset DISPLAY
#java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar LearnModel \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/binary \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/state_2 \
#  2 \
#  hg38 


### Run Forward and reverse separately
# Forward
#java -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar BinarizeBam \
#  /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/state_2_forward/ttseq_2state_forward_bam.txt \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/binary_forward

#unset DISPLAY
#java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar LearnModel \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/binary_forward \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/state_2_forward \
#  2 \
#  hg38 

# Reverse
java -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar BinarizeBam \
  /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/state_2_reverse/ttseq_2state_reverse_bam.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/binary_reverse

unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar LearnModel \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/binary_reverse \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/state_2_reverse \
  2 \
  hg38
