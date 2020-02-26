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

######################################################################
######################################################################
######################################################################
######### Code ran
## Reorder based on hierarchial clustering and re-label
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar Reorder \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/ReOrder/stateorderingfile.txt \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/ReOrder/labelmappingfile.txt \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/jurkat_15_segments.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/ReOrder/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/model_15.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/ReOrder

## Add the "other data"
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/ReOrder/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/ReOrder/overlap_dbs_all

######################################################################
######################################################################
######################################################################
### Reorder after Ernst and Kellis 2017
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar Reorder \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/stateorderingfile.txt \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/jurkat_15_segments.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/model_15.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst

unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/overlap_dbs_all2




### Now that Jurkat data is relabeled following Ernst
### Make overlap plots (for supplemental) against E115, E116, and E123

# Fix names for cells E115, E116, and E123
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115/E115_1_TssA.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed/A_E115_Active_TSS.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115/E115_2_TssAFlnk.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed/B_E115_Flanking_active_TSS.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115/E115_3_TxFlnk.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed/C_E115_Transcription_TSS_TTS.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115/E115_4_Tx.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed/D_E115_Strong_transcription.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115/E115_5_TxWk.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed/E_E115_Weak_transcription.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115/E115_6_EnhG.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed/F_E115_Genic_enhancer.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115/E115_7_Enh.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed/G_E115_Active_enhancer.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115/E115_8_ZNF.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed/H_E115_ZNF_genes___repeats.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115/E115_9_Het.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed/I_E115_Heterochromatin.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115/E115_10_TssBiv.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed/J_E115_Bivalent_poised_TSS.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115/E115_11_BivFlnk.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed/K_E115_Flanking_bivalent_TSS_enhancer.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115/E115_12_EnhBiv.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed/L_E115_Bivalent_enhancer.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115/E115_13_ReprPC.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed/M_E115_Repressed_Polycomb.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115/E115_14_ReprPCWk.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed/N_E115_Weak_repressed_Polycomb.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115/E115_15_Quies.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed/O_E115_Quiescent_low.bed

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116/E116_1_TssA.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed/A_E116_Active_TSS.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116/E116_2_TssAFlnk.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed/B_E116_Flanking_active_TSS.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116/E116_3_TxFlnk.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed/C_E116_Transcription_TSS_TTS.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116/E116_4_Tx.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed/D_E116_Strong_transcription.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116/E116_5_TxWk.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed/E_E116_Weak_transcription.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116/E116_6_EnhG.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed/F_E116_Genic_enhancer.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116/E116_7_Enh.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed/G_E116_Active_enhancer.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116/E116_8_ZNF.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed/H_E116_ZNF_genes___repeats.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116/E116_9_Het.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed/I_E116_Heterochromatin.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116/E116_10_TssBiv.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed/J_E116_Bivalent_poised_TSS.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116/E116_11_BivFlnk.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed/K_E116_Flanking_bivalent_TSS_enhancer.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116/E116_12_EnhBiv.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed/L_E116_Bivalent_enhancer.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116/E116_13_ReprPC.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed/M_E116_Repressed_Polycomb.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116/E116_14_ReprPCWk.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed/N_E116_Weak_repressed_Polycomb.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116/E116_15_Quies.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed/O_E116_Quiescent_low.bed

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123/E123_1_TssA.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed/A_E123_Active_TSS.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123/E123_2_TssAFlnk.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed/B_E123_Flanking_active_TSS.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123/E123_3_TxFlnk.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed/C_E123_Transcription_TSS_TTS.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123/E123_4_Tx.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed/D_E123_Strong_transcription.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123/E123_5_TxWk.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed/E_E123_Weak_transcription.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123/E123_6_EnhG.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed/F_E123_Genic_enhancer.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123/E123_7_Enh.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed/G_E123_Active_enhancer.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123/E123_8_ZNF.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed/H_E123_ZNF_genes___repeats.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123/E123_9_Het.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed/I_E123_Heterochromatin.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123/E123_10_TssBiv.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed/J_E123_Bivalent_poised_TSS.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123/E123_11_BivFlnk.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed/K_E123_Flanking_bivalent_TSS_enhancer.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123/E123_12_EnhBiv.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed/L_E123_Bivalent_enhancer.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123/E123_13_ReprPC.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed/M_E123_Repressed_Polycomb.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123/E123_14_ReprPCWk.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed/N_E123_Weak_repressed_Polycomb.bed
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123/E123_15_Quies.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed/O_E123_Quiescent_low.bed







mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/Supplemental_OverlapEnrichment_E115
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115_fixed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/Supplemental_OverlapEnrichment_E115/overlap_E115

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/Supplemental_OverlapEnrichment_E116
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116_fixed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/Supplemental_OverlapEnrichment_E116/overlap_E116

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/Supplemental_OverlapEnrichment_E123
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123_fixed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/Supplemental_OverlapEnrichment_E123/overlap_E123




















