#! /bin/bash
#SBATCH --job-name=Erose
#SBATCH --partition=super
#SBATCH --output=Erose.%j.out
#SBATCH --error=Erose.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### This script uses a 4 state hidden markov model to find transcribed/not in all ttseq, fwd ttseq, rev ttseq
### Then remove protein coding genes plus/minus 2kb
### Remove all RNAseq transcribed ("stable") regions
### Overlap with ChiP-seq
### Then with H3K4me1, H3K27ac, and H3K4me3 call enhancers with rose

module load samtools
module load python/2.7.6-epd
module load bedtools/2.25.0

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
###################################################################### ### Make DataBases
########## TTseq
### Filter 4 state models for transcribed
awk '{OFS="\t"} NR>1 && $4!=1 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/state_4/jurkat_4_999_dense.bed | sort -k 1,1 -k 2,2n | grep -v "_" | bedtools merge -d 5 -i stdin >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_AllTranscribed.bed

awk '{OFS="\t"} NR>1 && $4!=4 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/state_4_forward/jurkat_4_999_dense.bed | sort -k 1,1 -k 2,2n | grep -v "_" | bedtools merge -d 5 -i stdin >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_ForwardTranscribed.bed

awk '{OFS="\t"} NR>1 && $4!=1{print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/state_4_reverse/jurkat_4_999_dense.bed | sort -k 1,1 -k 2,2n | grep -v "_" | bedtools merge -d 5 -i stdin >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_ReverseTranscribed.bed

### Find all transcribed that contains both fwd and reverse
bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_AllTranscribed.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_ForwardTranscribed.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_Alltranscribed_inForward.bed

bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_AllTranscribed.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_ReverseTranscribed.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_Alltranscribed_inReverse.bed

bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_Alltranscribed_inForward.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_Alltranscribed_inReverse.bed | sort -k 1,1 -k 2,2n | uniq | grep -v "_" >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_AllTranscribed_inFandR.bed 



### Find F and R near each other but not in all
bedtools intersect -v -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_ForwardTranscribed.bed  -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_AllTranscribed.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/not_forward.bed

bedtools intersect -v -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_ReverseTranscribed.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_AllTranscribed.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/not_reverse.bed

cat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/not_forward.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/not_reverse.bed | sort -k 1,1 -k 2,2n | mergeBed -i stdin >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/not_merged.bed

bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/not_merged.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/not_forward.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/not_intersectA.bed
bedtools intersect -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/not_intersectA.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/not_reverse.bed | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_notAll_has_FandR.bed

rm not_forward.bed not_reverse.bed not_merged.bed not_intersectA.bed


### Find A,F,and R not in above 2:
cat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_AllTranscribed.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_ForwardTranscribed.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_ReverseTranscribed.bed | sort -k 1,1 -k 2,2n | bedtools merge -i stdin >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_AllFandR.bed

bedtools intersect -v -a <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_AllFandR.bed) -b <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_AllTranscribed_inFandR.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/delete.bed

bedtools intersect -v -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/delete.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_notAll_has_FandR.bed > /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_InOnly1_AllFandR.bed

rm /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/delete.bed




########## Get DB of protein coding genes +/- 2kb
awk '{OFS="\t"} $3=="gene" && $12 ~ /"protein_coding"/  && $7=="+" {print $1,$4-2000,$5+2000,"protein_coding",".",$7}' /project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/pc_plus.bed
awk '{OFS="\t"} $3=="gene" && $12 ~ /"protein_coding"/  && $7=="-" {print $1,$4-2000,$5+2000,"protein_coding",".",$7}' /project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/pc_minus.bed

cat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/pc_plus.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/pc_minus.bed | sort -k 1,1 -k 2,2n | grep "^chr" >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/protein_coding_PM2kb_toremove.bed

rm /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/pc_plus.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/pc_minus.bed


########## Get other annotated regions from RNAseq; to be removed
### The fpkm for each rep must be >1
### We will use bedtools subtract here
awk '{OFS="\t"} NR>1 {if ($4> 1 && $5> 1 && $6> 1  && $3 != "protein_coding") {print $1}}' /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq_reg/workflow/output_PE/countTable.fpkm.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/delete_db2remove.txt

grep -Ff /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/delete_db2remove.txt /project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf | awk '{OFS="\t"} NR>5 {if ($3 == "gene") {print $1, $4, $5, $10, $7, $8}}' >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/other_annotated_regions_toremove.bed

rm /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/delete_db2remove.txt


########## Remove the above databases from TTseq
cat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/protein_coding_PM2kb_toremove.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/other_annotated_regions_toremove.bed | sort -k 1,1 -k 2,2n | bedtools merge -i stdin > /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/all_regions2remove.bed


bedtools intersect -v -a <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_AllTranscribed_inFandR.bed ) -b <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/all_regions2remove.bed) > /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_AllTranscribed_inFandR_NotProteinCoding_NotTranscribedAnnotated.bed


bedtools intersect -v -a <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_notAll_has_FandR.bed) -b <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/all_regions2remove.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_notAll_has_FandR_NotProteinCoding_NotTranscribedAnnotated.bed


bedtools intersect -v -a <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_InOnly1_AllFandR.bed) -b <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/all_regions2remove.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_InOnly1_AllFandR_NotProteinCoding_NotTranscribedAnnotated.bed


######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
###################################################################### ### Make Super Enhancers
########## 1. Run rose - Combine before superEnhancers
cat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_AllTranscribed_inFandR_NotProteinCoding_NotTranscribedAnnotated.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_notAll_has_FandR_NotProteinCoding_NotTranscribedAnnotated.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_InOnly1_AllFandR_NotProteinCoding_NotTranscribedAnnotated.bed | awk '{OFS="\t"} {print $1, "ttseq_"NR, ".", $2, $3, ".", ".", ".", "ttseq_"NR}' >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/all_ttseq_possible_enhancers.gtf

# H3K27ac
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/H3K27ac

python /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/rose/rose/ROSE_main.py \
  -g HG18 \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/all_ttseq_possible_enhancers.gtf \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_merged.bam \
  -c /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Input_H3K27ac_merged.bam \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/H3K27ac \
  -s 12500 \
  -t 2500

### H3K4me3
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/H3K4me3

python /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/rose/rose/ROSE_main.py \
  -g HG18 \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/all_ttseq_possible_enhancers.gtf \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/H3K4me3_merge.bam \
  -c /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/Input_H3K4me3_merged.bam \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/H3K4me3 \
  -s 12500 \
  -t 2500


### H3K4me1
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/H3K4me1

python /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/rose/rose/ROSE_main.py \
  -g HG18 \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/all_ttseq_possible_enhancers.gtf \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/H3K4me1_single.bam \
  -c /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/Input_H3K4me1_single.bam \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/H3K4me1 \
  -s 12500 \
  -t 2500



### ttseq vs RNAseq
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ttseq_vs_rnaseq

python /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/rose/rose/ROSE_main.py \
  -g HG18 \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/all_ttseq_possible_enhancers.gtf \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Ttseq_CRAMER_merged.bam \
  -c /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/RNAseq.bam \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ttseq_vs_rnaseq \
  -s 12500 \
  -t 2500

### Filter out everything from gtf from ttseq superenhancers; ttseq was really good at calling links
awk '{OFS="\t"} NR>5 {if ($3 == "gene") {print $1, $4, $5, $10, $7, $8}}' /project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ttseq_vs_rnaseq/gencode_filter.bed

bedtools intersect -v -a <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ttseq_vs_rnaseq/all_ttseq_possible_enhancers_Gateway_SuperEnhancers.bed) -b <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ttseq_vs_rnaseq/gencode_filter.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ttseq_vs_rnaseq/ttseq_filteredSE.bed

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
###################################################################### ### Make Just Enhancers
########## Overlap each ttseq input (3) minus superEnhancer region with 3 histones

# SE; 767 total
cat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/ZP_H3K27ac_SE.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/ZQ_H3K4me3_SE.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/ZS_H3K4me1_SE.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/ZT_TTseq_SE.bed | sort -k 1,1 -k 2,2n | bedtools merge -i stdin >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/SuperEnhancers.bed

# remove SE from ttseq; 18357 left
bedtools intersect -v -a <(cat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_AllTranscribed_inFandR_NotProteinCoding_NotTranscribedAnnotated.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_notAll_has_FandR_NotProteinCoding_NotTranscribedAnnotated.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_InOnly1_AllFandR_NotProteinCoding_NotTranscribedAnnotated.bed | sort -k 1,1 -k 2,2n) -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/SuperEnhancers.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/ttseq_possible_enhancers_minusSuperEnhancers.bed

### Overlap with 3 histones
bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/ttseq_possible_enhancers_minusSuperEnhancers.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/consensusPeaks/H3K27ac.replicated.narrowPeak | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/H3K27ac_enhancers.bed

bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/ttseq_possible_enhancers_minusSuperEnhancers.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/consensusPeaks/H3K4me3.replicated.narrowPeak | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/H3K4me3_enhancers.bed

bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/ttseq_possible_enhancers_minusSuperEnhancers.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/GSM1603225_GFP_H3K4me1.replicated.narrowPeak | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/H3K4me1_enhancers.bed

# both f and rev
bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/ttseq_possible_enhancers_minusSuperEnhancers.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_ForwardTranscribed.bed | bedtools intersect -wa -a stdin -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/mkdb/jurkat_4_ReverseTranscribed.bed | sort -k 1,1 -k 2,2n | bedtools merge -i stdin | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/TTseq_enhancers.bed

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
###################################################################### ### Make Final Datasets
########## chromHMM and metagene plots; This is overlapping the 15 state chromHMM with the enhancer databases created
cat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/H3K27ac_enhancers.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/H3K4me3_enhancers.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/H3K4me1_enhancers.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/enhancers/TTseq_enhancers.bed | sort -k 1,1 -k 2,2n | bedtools merge -i stdin >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_enhancers.bed

cat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/H3K27ac/all_ttseq_possible_enhancers_Gateway_SuperEnhancers.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/H3K4me3/all_ttseq_possible_enhancers_Gateway_SuperEnhancers.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/H3K4me1/all_ttseq_possible_enhancers_Gateway_SuperEnhancers.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ttseq_vs_rnaseq/ttseq_filteredSE.bed | sort -k 1,1 -k 2,2n | bedtools merge -i stdin > /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed

# Overlap plots
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/overlap_refseq_db

# Neighborhood plots
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar NeighborhoodEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /work/BICF/s185797/programs/ChromHMM/ANCHORFILES/hg38/RefSeqTSS.hg38.txt.gz \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/neighborhood_refseqTSS

unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar NeighborhoodEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /work/BICF/s185797/programs/ChromHMM/ANCHORFILES/hg38/RefSeqTES.hg38.txt.gz \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/neighborhood_refseqTTS


########## Make metagene plots
module load deeptools/2.5.0.1
########## Enhancers
### ttseq
### merge fwd and rev bam files, and all bam files, and Convert to bigwig (in other script)
computeMatrix reference-point -S /archive/shared/DOrso_BICF/WashU_browser_data/ttseq_merged.bw /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_fwd.bw /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_rev.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_Enhancer.bed \
  --referencePoint center \
  -a 3000 \
  -b 3000 \
  -p max/2 \
  --sortRegions descend \
  --sortUsing mean \
  --sortUsingSamples 1 \
  --outFileSortedRegions ttseq_enhancers_heatmaps_order.bed \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/TTseq_enhancers.mat.gz

plotHeatmap -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/TTseq_enhancers.mat.gz \
  --samplesLabel TTseq_all TTseq_fwd TTseq_rev \
  --colorList 'white,red' \
  --refPointLabel center \
  --regionsLabel enhancers \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/TTseq_enhancers.pdf


### H3K27ac
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K27ac_pooled.fc_signal.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_Enhancer.bed \
  --referencePoint center \
  -a 3000 \
  -b 3000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K27_enhancers.mat.gz

plotHeatmap -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K27_enhancers.mat.gz \
  --samplesLabel H3K27ac \
  --colorList 'white,black' \
  --refPointLabel center \
  --regionsLabel enhancers \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K27ac_enhancers.pdf


### H3K4me3
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K4me3_pooled.fc_signal.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_Enhancer.bed \
  --referencePoint center \
  -a 3000 \
  -b 3000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K4me3_enhancers.mat.gz

plotHeatmap -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K27_enhancers.mat.gz \
  --samplesLabel H3K4me3 \
  --colorList 'white,black' \
  --refPointLabel center \
  --regionsLabel enhancers \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K4me3_enhancers.pdf

### H3K4me1
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/GSM1603225_GFP_H3K4me1_pooled.fc_signal.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_Enhancer.bed \
  --referencePoint center \
  -a 3000 \
  -b 3000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K4me1_enhancers.mat.gz

plotHeatmap -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K4me1_enhancers.mat.gz \
  --samplesLabel H3K4me1 \
  --colorList 'white,black' \
  --refPointLabel center \
  --regionsLabel enhancers \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K4me1_enhancers.pdf

### DNase
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/callPeaksMACS/DNAse-seq_ENCODE_GSM736501_pooled.fc_signal.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zr_Enhancer.bed \
  --referencePoint center \
  -a 3000 \
  -b 3000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/DNase_enhancers.mat.gz

plotHeatmap -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/DNase_enhancers.mat.gz \
  --samplesLabel DNase \
  --colorList 'white,blue' \
  --refPointLabel center \
  --regionsLabel enhancers \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/DNase_enhancers.pdf


#############
###############
#################
#
##################
########## Make metagene plots

########## SuperEnhancers
### ttseq
### merge fwd and rev bam files, and all bam files, and Convert to bigwig (in other script)
computeMatrix reference-point -S /archive/shared/DOrso_BICF/WashU_browser_data/ttseq_merged.bw /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_fwd.bw /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_rev.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed \
  --referencePoint center \
  -a 20000 \
  -b 20000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/TTseq_SuperEnhancers.mat.gz

plotHeatmap -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/TTseq_SuperEnhancers.mat.gz \
  --samplesLabel TTseq_all TTseq_fwd TTseq_rev \
  --colorList 'white,red' \
  --refPointLabel center \
  --regionsLabel SuperEnhancers \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/TTseq_SuperEnhancers.pdf


### H3K27ac
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K27ac_pooled.fc_signal.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed \
  --referencePoint center \
  -a 20000 \
  -b 20000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K27_SuperEnhancers.mat.gz

plotHeatmap -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K27_SuperEnhancers.mat.gz \
  --samplesLabel H3K27ac \
  --colorList 'white,black' \
  --refPointLabel center \
  --regionsLabel SuperEnhancers \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K27ac_SuperEnhancers.pdf


### H3K4me3
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K4me3_pooled.fc_signal.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed \
  --referencePoint center \
  -a 20000 \
  -b 20000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K4me3_SuperEnhancers.mat.gz

plotHeatmap -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K27_SuperEnhancers.mat.gz \
  --samplesLabel H3K4me3 \
  --colorList 'white,black' \
  --refPointLabel center \
  --regionsLabel SuperEnhancers \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K4me3_SuperEnhancers.pdf

### H3K4me1
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/GSM1603225_GFP_H3K4me1_pooled.fc_signal.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed \
  --referencePoint center \
  -a 20000 \
  -b 20000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K4me1_SuperEnhancers.mat.gz

plotHeatmap -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K4me1_SuperEnhancers.mat.gz \
  --samplesLabel H3K4me1 \
  --colorList 'white,black' \
  --refPointLabel center \
  --regionsLabel SuperEnhancers \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/H3K4me1_SuperEnhancers.pdf

### DNase
computeMatrix reference-point -S /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/callPeaksMACS/DNAse-seq_ENCODE_GSM736501_pooled.fc_signal.bw \
  -R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/zv_SuperEnhancer.bed \
  --referencePoint center \
  -a 20000 \
  -b 20000 \
  -p max/2 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/DNase_SuperEnhancers.mat.gz

plotHeatmap -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/DNase_SuperEnhancers.mat.gz \
  --samplesLabel DNase \
  --colorList 'white,blue' \
  --refPointLabel center \
  --regionsLabel SuperEnhancers \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/four_state_chromHMM_TTseq/ttseq_rose_enhancers/ppdb/heatmaps/DNase_SuperEnhancers.pdf

