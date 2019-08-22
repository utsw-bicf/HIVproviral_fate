#! /bin/bash
#SBATCH --job-name=rose
#SBATCH --partition=super
#SBATCH --output=rose_rose.%j.out
#SBATCH --error=rose_rose.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### This script uses a 2 state hidden markov model to find transcribed/not in all ttseq, fwd ttseq, rev ttseq
### Then with H3K4me1, H3K27ac, and H3K4me3 call super enhancers with rose (Don't remove anything!)

module load samtools
module load python/2.7.6-epd
module load bedtools/2.25.0

########## TTseq
### Filter 2 state models for transcribed
awk '{OFS="\t"} NR>1 && $4==2 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/state_2/jurkat_2_dense.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_2_AllTranscribed.bed

awk '{OFS="\t"} NR>1 && $4==2 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/state_2_forward/jurkat_2_dense.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_2_ForwardTranscribed.bed

awk '{OFS="\t"} NR>1 && $4==2 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/two_state_chromHMM_TTseq/state_2_reverse/jurkat_2_dense.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_2_ReverseTranscribed.bed

### Find all transcribed that contains both fwd and reverse
bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_2_AllTranscribed.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_2_ForwardTranscribed.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_2_Alltranscribed_inForward.bed

bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_2_AllTranscribed.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_2_ReverseTranscribed.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_2_Alltranscribed_inReverse.bed

bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_2_Alltranscribed_inForward.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_2_Alltranscribed_inReverse.bed | sort -k 1,1 -k 2,2n | uniq | grep -v "_" >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_AllTranscribed_inFandR.bed 

# Make gff file for rose
awk '{OFS="\t"} {print $1, "ttseq_"NR, ".", $2, $3, ".", ".", ".", "ttseq_"NR}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_AllTranscribed_inFandR.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_AllTranscribed_inFandR.gff


########## Make histone databases
### H3K27ac
#samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_merged.bam ${loc}/SRR2043614.filt.nodup.bam ${loc}/SRR1603650.filt.nodup.bam ${loc}/SRR1603654.filt.nodup.bam ${loc}/GSM1603211.filt.nodup.bam
#samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_merged.bam

#samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Input_H3K27ac_merged.bam ${loc}/SRR2043612.filt.nodup.bam ${loc}/SRR1603649.filt.nodup.bam ${loc}/SRR1603652.filt.nodup.bam ${loc}/GSM1603229.filt.nodup.bam
#samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Input_H3K27ac_merged.bam

#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Input_H3K27ac_merged.bam .
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Input_H3K27ac_merged.bam.bai

### H3K4me1
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603225.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/H3K4me1_single.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603225.filt.nodup.bam.bai /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/H3K4me1_single.bam.bai

#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603229.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/Input_H3K4me1_single.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603229.filt.nodup.bam.bai /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/Input_H3K4me1_single.bam.bai

### H3K4me3
#samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K4me3_consensus.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/filterReads/SRR577482.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/filterReads/SRR577483.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/filterReads/GSM1603213.filt.nodup.bam
#samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K4me3_consensus.bam

#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/H3K4me3_consensus.bam /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/H3K4me3_merge.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/H3K4me3_consensus.bam.bai /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/H3K4me3_merge.bam.bai

#samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/Input_H3K4me3_merged.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/SRR577484.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603229.filt.nodup.bam
#samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/Input_H3K4me3_merged.bam

########## Run Rose
### Must be run from Rose folder, which is at
### /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/rose/rose
python /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/rose/rose/ROSE_main.py \
  -g HG18 \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_AllTranscribed_inFandR.gff \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_merged.bam \
  -c /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Input_H3K27ac_merged.bam \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_SEresults \
  -s 12500 \
  -t 2500


python /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/rose/rose/ROSE_main.py \
  -g HG18 \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_AllTranscribed_inFandR.gff \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/H3K4me1_single.bam \
  -c /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/Input_H3K4me1_single.bam \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_SEresults \
  -s 12500 \
  -t 2500


python /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/rose/rose/ROSE_main.py \
  -g HG18 \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_AllTranscribed_inFandR.gff \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/H3K4me3_merge.bam \
  -c /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/Input_H3K4me3_merged.bam \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_SEresults \
  -s 12500 \
  -t 2500



########## Make enhancer and super enhancer db
### H3K4me1
awk '{OFS="\t"} $10 == 1 {print $2,$3,$4,$1,$9,"."}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_results/jurkat_AllTranscribed_inFandR_AllEnhancers.table.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_results/H3K4me1_SE.bed

cp /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_results/H3K4me1_SE.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/

### H3K27ac
awk '{OFS="\t"} $10 == 1 {print $2,$3,$4,$1,$9,"."}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_results/jurkat_AllTranscribed_inFandR_AllEnhancers.table.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_results/H3K27ac_SE.bed

cp /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_results/H3K27ac_SE.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/

### H3K4me3
awk '{OFS="\t"} $10 == 1 {print $2,$3,$4,$1,$9,"."}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_results/jurkat_AllTranscribed_inFandR_AllEnhancers.table.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_results/H3K4me3_SE.bed

cp /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_results/H3K4me3_SE.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/


########## Run intervene to do upset plot
intervene upset -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/H3K27ac_SE.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/H3K4me1_SE.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/H3K4me3_SE.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/dbSuper_Jurkat_hg19_LO_hg38.bed\
  --type genomic \
  --names H3K27ac,H3K4me1,H3K4me3,dbSuper \
  --project SuperEnhancers_upset

Rscript /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/chromHMM_results/Intervene_results/SuperEnhancers_upset_upset.R



