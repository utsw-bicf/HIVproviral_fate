#! /bin/bash
#SBATCH --job-name=Erose
#SBATCH --partition=super
#SBATCH --output=Erose.%j.out
#SBATCH --error=Erose.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### This script uses a 2 state hidden markov model to find transcribed/not in all ttseq, fwd ttseq, rev ttseq
### Then remove protein coding genes plus/minus 2kb
### Remove all RNAseq transcribed ("stable") regions
### Overlap with ChiP-seq
### Then with H3K4me1, H3K27ac, and H3K4me3 call enhancers with rose

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



########## Get DB of protein coding genes +/- 2kb
awk '{OFS="\t"} $3=="gene" && $12 ~ /"protein_coding"/  && $7=="+" {print $1,$4-2000,$5+2000,"protein_coding",".",$7}' /project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/pc_plus.bed
awk '{OFS="\t"} $3=="gene" && $12 ~ /"protein_coding"/  && $7=="-" {print $1,$4-2000,$5+2000,"protein_coding",".",$7}' /project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/pc_minus.bed

cat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/pc_plus.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/pc_minus.bed | sort -k 1,1 -k 2,2n | grep "^chr" >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/protein_coding_PM2kb_toremove.bed

rm /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/pc_plus.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/pc_minus.bed


########## Get other annotated regions from RNAseq; to be removed
### The fpkm for each rep must be >1
### We will use bedtools subtract here
awk '{OFS="\t"} NR>1 {if ($4> 1 && $5> 1 && $6> 1  && $3 != "protein_coding") {print $1}}' /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq_reg/workflow/output_PE/countTable.fpkm.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/delete_db2remove.txt

grep -Ff /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/delete_db2remove.txt /project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf | awk '{OFS="\t"} NR>5 {if ($3 == "gene") {print $1, $4, $5, $10, $7, $8}}' >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/other_annotated_regions_toremove.bed

rm /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/delete_db2remove.txt


########## Remove the above databases from TTseq
bedtools subtract -a <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_AllTranscribed_inFandR.bed) -b <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/protein_coding_PM2kb_toremove.bed) > /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_AllTranscribed_inFandR_NotProteinCoding.bed

bedtools subtract -a <(sort -k 1,1 -k 2,2n  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_AllTranscribed_inFandR_NotProteinCoding.bed) -b <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/other_annotated_regions_toremove.bed) > /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_AllTranscribed_inFandR_NotProteinCoding_NotTranscribedAnnotated.bed



########## 1. Run rose
awk '{OFS="\t"} {print $1, "ttseq_"NR, ".", $2, $3, ".", ".", ".", "ttseq_"NR}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_AllTranscribed_inFandR_NotProteinCoding_NotTranscribedAnnotated.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/rose_Number1.gff

python /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/rose/rose/ROSE_main.py \
  -g HG18 \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/rose_Number1.gff \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_merged.bam \
  -c /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Input_H3K27ac_merged.bam \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/Number1 \
  -s 12500 \
  -t 2500

awk '{OFS="\t"} $10 == 0 {print $2,$3,$4,$1,$9,"."}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/Number1/rose_Number1_AllEnhancers.table.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/Number1/number1_enhancers.bed



########## 2. Find chip peaks that overlap with ttseq then run rose
bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_AllTranscribed_inFandR_NotProteinCoding_NotTranscribedAnnotated.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/consensusPeaks/H3K27ac.replicated.narrowPeak >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_AllTranscribed_inFandR_NotProteinCoding_NotTranscribedAnnotated_hasH3K27ac.bed

awk '{OFS="\t"} {print $1, "ttseq_"NR, ".", $2, $3, ".", ".", ".", "ttseq_"NR}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/jurkat_AllTranscribed_inFandR_NotProteinCoding_NotTranscribedAnnotated_hasH3K27ac.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/rose_Number2.gff

python /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/rose/rose/ROSE_main.py \
  -g HG18 \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/rose_Number2.gff \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_merged.bam \
  -c /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Input_H3K27ac_merged.bam \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/Number2 \
  -s 12500 \
  -t 2500

awk '{OFS="\t"} $10 == 0 {print $2,$3,$4,$1,$9,"."}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/Number2/rose_Number2_AllEnhancers.table.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/Number2/number2_enhancers.bed



########## 3. Filter previous run
awk '{OFS="\t"} $10 == 0 {print $2,$3,$4,$1,$9,"."}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_SEresults/jurkat_AllTranscribed_inFandR_AllEnhancers.table.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/Number3/H3K27ac_possible_enhancers.bed

### Remove protein coding
bedtools intersect -v -a <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/Number3/H3K27ac_possible_enhancers.bed) -b <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/protein_coding_PM2kb_toremove.bed) > /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/Number3/number3_nogenes.bed

### Remove "other" annotations
bedtools intersect -v -a <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/Number3/number3_nogenes.bed) -b <(sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/other_annotated_regions_toremove.bed) >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/Number3/number3_enhancers.bed




########## 4. Run overlap plot
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/chromHMM_results/enhancers_3ways


########## Run intervene to do upset plot
intervene upset -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/number*.bed\
  --type genomic \
  --names number1,number2,number3 \
  --project enhancers_3ways

module load R
Rscript /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/chromHMM_results/Intervene_results/enhancers_3ways_upset.R


































########## Filter rose output
### remove overlapping genes regions + 2kb 
awk '{OFS="\t"} $3=="gene" && $12 ~ /"protein_coding"/  && $7=="+" {print $1,$4-2000,$5,"protein_coding",".",$7}' /project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/pc_plus.bed
awk '{OFS="\t"} $3=="gene" && $12 ~ /"protein_coding"/  && $7=="-" {print $1,$4,$5+2000,"protein_coding",".",$7}' /project/shared/bicf_workflow_ref/human/GRCh38/gencode.gtf >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/pc_minus.bed

cat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/pc_plus.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/pc_minus.bed | sort -k 1,1 -k 2,2n | grep "^chr" >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/genes_and_promoters.bed

rm /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/pc_plus.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/pc_minus.bed



########## Make enhancer and super enhancer db
### H3K4me1
awk '{OFS="\t"} $10 == 1 {print $2,$3,$4,$1,$9,"."}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_results/jurkat_AllTranscribed_inFandR_AllEnhancers.table.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_results/H3K4me1_SE.bed

awk '{OFS="\t"} $10 == 0 {print $2,$3,$4,$1,$9,"."}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_results/jurkat_AllTranscribed_inFandR_AllEnhancers.table.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_results/H3K4me1_possible_enhancers.bed

bedtools intersect -wo -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_results/H3K4me1_possible_enhancers.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/genes_and_promoters.bed | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_results/H3K4me1_possible_enhancers_annotated.bed

awk '{OFS="\t"} $13/($9-$8) >= 0.5 {print $1,$2,$3,$4,$5,$6}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_results/H3K4me1_possible_enhancers_annotated.bed | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_results/enhancers_in_genes_and_promoters.bed

bedtools intersect -v -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_results/H3K4me1_possible_enhancers.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_results/enhancers_in_genes_and_promoters.bed | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me1_Rose_results/H3K4me1_enhancers.bed

### H3K27ac
awk '{OFS="\t"} $10 == 1 {print $2,$3,$4,$1,$9,"."}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_results/jurkat_AllTranscribed_inFandR_AllEnhancers.table.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_results/H3K27ac_SE.bed

awk '{OFS="\t"} $10 == 0 {print $2,$3,$4,$1,$9,"."}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_results/jurkat_AllTranscribed_inFandR_AllEnhancers.table.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_results/H3K27ac_possible_enhancers.bed

bedtools intersect -wo -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_results/H3K27ac_possible_enhancers.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/genes_and_promoters.bed | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_results/H3K27ac_possible_enhancers_annotated.bed

awk '{OFS="\t"} $13/($9-$8) >= 0.5 {print $1,$2,$3,$4,$5,$6}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_results/H3K27ac_possible_enhancers_annotated.bed | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_results/enhancers_in_genes_and_promoters.bed

bedtools intersect -v -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_results/H3K27ac_possible_enhancers.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_results/enhancers_in_genes_and_promoters.bed | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K27ac_Rose_results/H3K27ac_enhancers.bed


### H3K4me3
awk '{OFS="\t"} $10 == 1 {print $2,$3,$4,$1,$9,"."}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_results/jurkat_AllTranscribed_inFandR_AllEnhancers.table.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_results/H3K4me3_SE.bed

awk '{OFS="\t"} $10 == 0 {print $2,$3,$4,$1,$9,"."}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_results/jurkat_AllTranscribed_inFandR_AllEnhancers.table.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_results/H3K4me3_possible_enhancers.bed

bedtools intersect -wo -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_results/H3K4me3_possible_enhancers.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/makeDB/genes_and_promoters.bed | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_results/H3K4me3_possible_enhancers_annotated.bed

awk '{OFS="\t"} $13/($9-$8) >= 0.5 {print $1,$2,$3,$4,$5,$6}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_results/H3K4me3_possible_enhancers_annotated.bed | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_results/enhancers_in_genes_and_promoters.bed

bedtools intersect -v -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_results/H3K4me3_possible_enhancers.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_results/enhancers_in_genes_and_promoters.bed | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/ttseq_H3K4me3_Rose_results/H3K4me3_enhancers.bed





########## Run chromHMM

unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/chromHMM_results/new_rose_enhancers


########## Run intervene to do upset plot
intervene upset -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/H3K27ac_SE.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/H3K4me1_SE.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/H3K4me3_SE.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/dbSuper_Jurkat_hg19_LO_hg38.bed\
  --type genomic \
  --names H3K27ac,H3K4me1,H3K4me3,dbSuper \
  --project SuperEnhancers_upset

Rscript /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/chromHMM_results/Intervene_results/SuperEnhancers_upset_upset.R


intervene upset -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/H3K27ac_enhancers.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/H3K4me1_enhancers.bed /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/refseq_annotations/H3K4me3_enhancers.bed \
  --type genomic \
  --names H3K27ac,H3K4me1,H3K4me3 \
  --project Enhancers_upset

Rscript /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ttseq_histone_Rose/chromHMM_results/Intervene_results/Enhancers_upset_upset.R
