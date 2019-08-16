#!/bin/bash

#SBATCH --job-name=TTeDB
#SBATCH --partition=super
#SBATCH --output=TTeDB.%j.out
#SBATCH --error=TTeDB.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

########## This script makes and enhancer/super enhancer database from TT seq data

module load bedtools
module load deeptools/2.5.0.1 
module load samtools
module load pybedtools
module load python/2.7.x-anaconda

### use plotProfile to find the cutoffs around the TSS and TES for promoter and polII falls off
computeMatrix scale-regions -S /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.fwd.bigWig /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.fwd.bigWig /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.rev.bw /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.rev.bw \
  -R /project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf \
  --beforeRegionStartLength 3000 \
  --regionBodyLength 5000 \
  --afterRegionStartLength 3000 \
  --skipZeros \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/matrix.mat.gz

plotProfile -m matrix.mat.gz -out TTseq_profile.pdf --numPlotsPerRow 1

########## Note the above graph isn't clear; need to instead of removing genes + 2Kb, remove genes +2Kb complete overlap, plus remove annotated lincRNAs

### Remove un-transcribed regions (1) from all TTseq chromHMM 2 state
awk '{OFS="\t"} NR>1 && $4==2 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/state_2/jurkat_2_dense.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_transcribed.bed

### Remove annotated features plus 2kb upstream of TSS from transcribed regions
# Fix gtf file
# protein-coding gene, add promoter
awk '{OFS="\t"} $3=="gene" && $12 ~ /"protein_coding"/  && $7=="+" {print $1,$4-2000,$5,"protein_coding",".",$7}' /project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/pc_plus.bed
awk '{OFS="\t"} $3=="gene" && $12 ~ /"protein_coding"/  && $7=="-" {print $1,$4,$5+2000,"protein_coding",".",$7}' /project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/pc_minus.bed
grep 'transcript_type \"lincRNA\"' /project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf | awk '{OFS="\t"} {print $1, $4, $5, $3, ".", $7}' >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/linc.bed 

cat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/pc_plus.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/pc_minus.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/linc.bed | sort -k 1,1 -k 2,2n | grep "^chr" >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/known_transcript_db2remove.bed


rm /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/pc_plus.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/pc_minus.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/linc.bed

bedtools intersect -v -a <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_transcribed.bed) -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/known_transcript_db2remove.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed.bed

#bedtools subtract -a <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_transcribed.bed) -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/known_transcript_db2remove.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed.bed

### Find regions that contain both the + and - TTseq strands
awk '{OFS="\t"} NR>1 && $4==2 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/state_2_forward/jurkat_2_dense.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_ForwardTranscribed.bed
awk '{OFS="\t"} NR>1 && $4==2 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/state_2_reverse/jurkat_2_dense.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_ReverseTranscribed.bed

bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_ForwardTranscribed.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inForward.bed
bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_ReverseTranscribed.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inReverse.bed

bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inForward.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inReverse.bed | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR.bed 

### Annotate each region, then annotate by chromHMM group
bedtools intersect -wao -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR.bed \
-b <(awk 'NR>1 {print $0}' /project/shared/bicf_workflow_ref/human/GRCh38/genenames.txt) \
>/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_DELETE.bed

bedtools intersect -wao -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_DELETE.bed \
-b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed | awk '{OFS="\t"} $1 !~/_/ {print $0}' >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_annotated.bed

rm /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_DELETE.bed

# Split into able to annotate and not
awk '{OFS="\t"} $11 ~ /-1/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_annotated.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_partially_annotated.bed
awk '{OFS="\t"} $11 !~ /-1/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_annotated.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_fully_annotated.bed















### Merge "enhancers" that are within 12500 of each other
bedtools merge -d 12500 -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR.bed | grep -v "_" >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_merged12500.bed

### calculate rpkm
echo "experiment,file
ttfwd87,/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.fwd.bam
ttfwd88,/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.fwd.bam
ttrev87,/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.rev.bam
ttrev88,/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.rev.bam" >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/calc_rpkm.csv

python /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/machine_learning/venkat_rpkm_3.py \
--peaks /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_merged12500.bed \
--experiments /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/calc_rpkm.csv \
-f jurkat_2_notgenes_transcribed_inFandR_merged12500_rpkm \
--minimum 0
 
### Filter enhancers for > 150 bp, and then total rpkm >2,10,50 and run overlap 
# 150 only
awk '{OFS="\t"} NR>1 {if ($3-$2+1 > 150) print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_merged12500_rpkm_filtered_peaks.tsv >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/enhancers_150bp.bed
# 150 and >2
awk '{OFS="\t"} NR>1 {if ($3-$2+1 > 150 && $4+$5+$6+$7 > 2) print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_merged12500_rpkm_filtered_peaks.tsv >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/enhancers_150bp_2rpkm.bed
# 150 bp and >10
awk '{OFS="\t"} NR>1 {if ($3-$2+1 > 150 && $4+$5+$6+$7 > 10) print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_merged12500_rpkm_filtered_peaks.tsv >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/enhancers_150bp_10rpkm.bed
# 150 bp and >50
awk '{OFS="\t"} NR>1 {if ($3-$2+1 > 150 && $4+$5+$6+$7 > 50) print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_merged12500_rpkm_filtered_peaks.tsv >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/enhancers_150bp_50rpkm.bed

### Filter >150 & <9000, and >9000
awk '{OFS="\t"} NR>1 {if ($3-$2+1 > 150 && $3-$2+1 < 9000) print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_merged12500_rpkm_filtered_peaks.tsv >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/enhancers_150to9000bp.bed

awk '{OFS="\t"} NR>1 {if ($3-$2+1 > 9000) print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_merged12500_rpkm_filtered_peaks.tsv >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/enhancers_gt9000bp.bed

# Run overlap plots to see if good
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/overlap_test





### Annotate "enhancers" 
# Make bed file of annotations
awk '{OFS="\t"}NR>5{print $1,$4,$5,$10,".",$7}' /project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/gencode.bed

bedtools intersect -wao -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR.bed \
-b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/gencode.bed \
>/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR_annotated.bed





################Dont use
########## Make sliding window of reference genome
### Make chrom sizes file
#bedtools makewindows \
#  -g /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes \
#  -w 12500 \
#  -s 200 \
#  >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/GRCh38_by12500_every200.bed


########## Make counts for possible enhancers
#bedtools intersect \
#  -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/GRCh38_by12500_every200.bed \
#  -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR.bed \
#  -c \
#  >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/make_db/jurkat_2_notgenes_transcribed_inFandR.countsin12500.bed



