#!/bin/bash

#SBATCH --job-name=chromHMM
#SBATCH --partition=super
#SBATCH --output=HMM.%j.out
#SBATCH --error=HMM.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

#### This script runs chromHMM histone plus MNase/DNAase and RNase data
#### For only 15 states

######################################################################
######################################################################
########## Start of Done previously
######################################################################
######################################################################

######### softlink all bam files to a single folder
### DNase-seq
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/filterReads/ENCFF001DPG.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/DNAse-seq_rep1.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/filterReads/ENCFF001DPF.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/DNAse-seq_rep2.bam
### MNase-seq
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/filterReads/Jkt106_MNase_Seq_20170218_1_val_1Jkt106_MNase_Seq_20170218_2_val_2.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/MNAse-seq.bam
### Chipseq
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR2043614.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/H3K27ac_YOUNG_GSM1697882.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR2043612.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/INPUT_YOUNG_GSM1697880.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR1603650.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/H3K27ac_YOUNG_GSM1519638.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR1603649.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/INPUT_YOUNG_GSM1519637.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR577482.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/H3K4me3_ENCODE_GSM945267_rep1.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR577483.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/H3K4me3_ENCODE_GSM945267_rep2.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR577484.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/INPUT_ENCODE_GSM945268.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR1603654.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/H3K27ac_YOUNG_GSM1519642.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR1603652.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/INPUT_YOUNG_GSM1519640.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR2157609.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/PolII_N20_YOUNG_GSM1850204.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR1057276.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/INPUT_YOUNG_GSM1296386.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR3722577.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/PolII_N20_YOUNG_GSM2218766.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR3722573.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/INPUT_YOUNG_GSM2218759.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603225.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/H3K4me1_DORSO_GSM1603225.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603229.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/INPUT_DORSO_GSM1603229.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603213.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/H3K4me3_DORSO_GSM1603213.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603209.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/H3K36me3_DORSO_GSM1603209.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603215.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/H3K79me3_DORSO_GSM1603215.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603215.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/H3K9me3_DORSO_GSM1603227.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603211.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/H3K27ac_DORSO_GSM1603211.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603217.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/PolII_4H8_DORSO_GSM1603217.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603223.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/PolII_S5P_DORSO_GSM1603223.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603221.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/PolII_S2P_DORSO_GSM1603221.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR3722566.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/BRD4_YOUNG_GSM2218755.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR3722570.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/INPUT_YOUNG_GGSM2218759.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR1522115.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/MED1_YOUNG_GSM1442004.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR1522118.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/INPUT_YOUNG_GSM1442007.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR6010201.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/YY1_YOUNG_GSM2773998.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR6010202.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/INPUT_YOUNG_GSM2773999.bam
### TT-seq
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/Ttseq_CRAMER_GSM2260187.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/Ttseq_CRAMER_GSM2260188.bam
### RNAseq
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq_reg/workflow/output_PE/GSM2260195.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/RNAseq_CRAMER_GSM2260195.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq_reg/workflow/output_PE/Emily_jurkat_A.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/Emily_jurkat_rep1.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq_reg/workflow/output_PE/Emily_jurkat_B.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/Emily_jurkat_rep2.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq_reg/workflow/output_PE/Emily_jurkat_C.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/Emily_jurkat_rep3.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq_reg/workflow/output_SE/SRR2130272.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/RNAseq_KEENE_GSM1833433.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq_reg/workflow/output_SE/SRR2130273.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/RNAseq_KEENE_GSM1833434.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq_reg/workflow/output_SE/SRR2130274.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/RNAseq_KEENE_GSM1833435.bam
### 4sU
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_4sU/S4sUseq_KEENE_GSM1833447.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/S4sUseq_KEENE_GSM1833447.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_4sU/S4sUseq_KEENE_GSM1833448.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/S4sUseq_KEENE_GSM1833448.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_4sU/S4sUseq_KEENE_GSM1833449.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/S4sUseq_KEENE_GSM1833449.bam
### CTCF
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/chipseq_analysis_4CTCF/workflow/output_CTCF/filterReads/SRR2029842.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/CTCF_YOUNG_GSM1689151.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/chipseq_analysis_4CTCF/workflow/output_CTCF/filterReads/SRR2029843.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/CTCF_YOUNG_GSM1689152.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/atacseq_analysis_4CTCFnoinput/workflow/output/filterReads/SRR014986.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/CTCF_ZHAO_GSM325899.bam

######################################################################
######################################################################
########## End of Done previously
######################################################################
######################################################################

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_and_ases/binary
java -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar BinarizeBam \
  /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/input_bam/histone_and_ases_input_bam.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_and_ases/binary

### Run the modeling
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar LearnModel \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_and_ases/binary \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_and_ases/state15 \
  15 \
  hg38 

### Annotate states
### Look in chromHMM_all.sh for annotation files
### Download files from https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/
### Ivan said that E115, E116, and E123 might be good, can also try all

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_and_ases/OverlapEnrichment_E115
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_and_ases/state15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E115 \
  OverlapEnrichment_E115/overlap_15_E115

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_and_ases/OverlapEnrichment_E116
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_and_ases/state15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E116 \
  OverlapEnrichment_E116/overlap_15_E116

mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_and_ases/OverlapEnrichment_E123
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_and_ases/state15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E123 \
  OverlapEnrichment_E117/overlap_15_E123


### Re-run with enhancer and superenhancer databases
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_and_ases/OverlapEnrichment_enhancers
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_and_ases/state15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_and_ases/OverlapEnrichment_enhancers/overlap_enhancers





















