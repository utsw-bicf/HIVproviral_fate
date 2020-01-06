#!/bin/bash

#SBATCH --job-name=chromHMM
#SBATCH --partition=super
#SBATCH --output=HMM.%j.out
#SBATCH --error=HMM.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

#### This script runs chromHMM histone plus MNase/DNAase, damIDseq and RNase data
#### For only 15 states

######### softlink all bam files to a single folder
### Everything but damIDseq previously done

### damID
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/counts/SRR5261759_filt.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/INPUT_DELASHERAS_SRR5261759.bam
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/counts/SRR5261760_filt.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/damID_DELASHERAS_SRR5261760.bam
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/counts/SRR5261761_filt.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/INPUT_DELASHERAS_SRR5261761.bam
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/counts/SRR5261762_filt.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/damID_DELASHERAS_SRR5261762.bam

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_ases_dam/binary
java -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar BinarizeBam \
  /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/histone_ases_dam_input_bam.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_ases_dam/binary

### Run the modeling
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar LearnModel \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_ases_dam/binary \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_ases_dam/state15 \
  15 \
  hg38 

### Annotate states
### Look in chromHMM_all.sh for annotation files
### Download files from https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/
### Ivan said that E115, E116, and E123 might be good, can also try all

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_ases_dam/OverlapEnrichment_E115
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_ases_dam/state15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115 \
  OverlapEnrichment_E115/overlap_15_E115

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_ases_dam/OverlapEnrichment_E116
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_ases_dam/state15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116 \
  OverlapEnrichment_E116/overlap_15_E116

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_ases_dam/OverlapEnrichment_E123
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_ases_dam/state15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123 \
  OverlapEnrichment_E117/overlap_15_E123


### Re-run with enhancer and superenhancer databases
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_ases_dam/OverlapEnrichment_enhancers
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_ases_dam/state15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_ases_dam/OverlapEnrichment_enhancers/overlap_enhancers






















