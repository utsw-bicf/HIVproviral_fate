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
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/binary
java -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar BinarizeBam \
  /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/histone_input_bam.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/binary

### Run the modeling
/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states
unset DISPLAY
for i in $(seq 5 25); do
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar LearnModel \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/binary \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn${i} \
  ${i} \
  hg38 
done

### Compare models
## softlink all but emissions_25 to a single folder
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/emissions
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn*/emissions_*.txt /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/emissions
rm /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/emissions/emissions_25.txt
java -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar CompareModels \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn25/emissions_25.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/emissions \
  histone_compare_model

### Annotate states
### Look in chromHMM_all.sh for annotation files
### Download files from https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/
### Ivan said that E115, E116, and E123 might be good, can also try all

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/OverlapEnrichment_E115
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115 \
  OverlapEnrichment_E115/overlap_15_E115

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/OverlapEnrichment_E116
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116 \
  OverlapEnrichment_E116/overlap_15_E116

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/OverlapEnrichment_E123
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123 \
  OverlapEnrichment_E117/overlap_15_E123

### To pick the correct number of states, run Rscript
### Or find a cell type that has states already done and overlap the results, with the OverlapEnrichment command
### Look at https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html
### or https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/

























