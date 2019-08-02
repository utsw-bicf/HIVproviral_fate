#!/bin/bash

#SBATCH --job-name=histone
#SBATCH --partition=super
#SBATCH --output=histone.%j.out
#SBATCH --error=histone.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

#### This script runs chromHMM on all data

######### softlink all bam files to a single folder
### Look in chromHMM_all.sh

### turn bams into binary bins
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/binary
java -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar BinarizeBam \
  /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/histone_input_bam7.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/binary

### Run the modeling
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar LearnModel \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/binary \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15 \
  15 \
  hg38 

## Add the "other data"
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/overlap_dbs

### Annotate states
### Look in chromHMM_all.sh for annotation files
### Download files from https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/
### Ivan said that E115, E116, and E123 might be good, can also try all

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/OverlapEnrichment_E115
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115 \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/OverlapEnrichment_E115/overlap_15_E115

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/OverlapEnrichment_E116
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116 \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/OverlapEnrichment_E116/overlap_15_E116

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/OverlapEnrichment_E123
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123 \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/OverlapEnrichment_E123/overlap_15_E123

### To pick the correct number of states, run Rscript
### Or find a cell type that has states already done and overlap the results, with the OverlapEnrichment command
### Look at https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html
### or https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/












######### Code ran
## Reorder and label
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar Reorder \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder/stateorderingfile.txt \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder/labelmappingfile.txt \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/jurkat_15_segments.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/model_15.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder

## Add the "other data"
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder/overlap_dbs

########## Ivan reordered 2
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar Reorder \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder2/stateorderingfile.txt \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder2/labelmappingfile.txt \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/jurkat_15_segments.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder2/jurkat_15_segments_reorder2.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/model_15.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder2

## Add the "other data"
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder2/jurkat_15_segments_reorder2.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder2/overlap_dbs


########## cluster reordered 3
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar Reorder \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder3_cluster/stateorderingfile.txt \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder3_cluster/labelmappingfile.txt \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/jurkat_15_segments.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder3_cluster/jurkat_15_segments_reorder3_cluster.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/model_15.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder3_cluster

## Add the "other data"
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder3_cluster/jurkat_15_segments_reorder3_cluster.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder3_cluster/overlap_dbs






















