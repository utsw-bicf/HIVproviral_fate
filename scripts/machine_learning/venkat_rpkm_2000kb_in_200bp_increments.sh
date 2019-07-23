#!/bin/bash

#SBATCH --job-name=HIV_rpkm
#SBATCH --partition=super
#SBATCH --output=HIV_rpkm.%j.out
#SBATCH --error=HIV_rpkm.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load pybedtools
module load python/2.7.x-anaconda
module load samtools
module load bedtools

### Make 200 bp bins HIV expression/insertion site
### Need to remove duplicates which are 1)chr3	48247567	48251569, 2)chr4	79730674	79734676
for i in {-10..9}; do
    left=$((${i}*200))
    right=$((${left}+200))
    awk -v left="$left" -v right="$right" '{OFS="\t"} {print $1, $2+left, $3+right, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.bed | grep -v "_" | sort -k 1,1 -k2,2n | uniq -u >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/hiv_expression_ML_${left}_${right}.bed
done

### Make csv of input bam files
### Merge bam files
### H3K27ac
#samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K27ac_consensus.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/filterReads/SRR2043614.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/filterReads/SRR1603650.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/filterReads/SRR1603654.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/filterReads/GSM1603211.filt.nodup.bam
#samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K27ac_consensus.bam

### H3K36me3
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603209.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K36me3_single.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603209.filt.nodup.bam.bai /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K36me3_single.bam.bai

### H3K79me3
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603215.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K79me3_single.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603215.filt.nodup.bam.bai /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K79me3_single.bam.bai

### H3K9me3
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603227.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K9me3_single.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603227.filt.nodup.bam.bai /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K9me3_single.bam.bai

### H3K4me1
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603225.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K4me1_single.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/filterReads/GSM1603225.filt.nodup.bam.bai /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K4me1_single.bam.bai

### H3K4me3
#samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K4me3_consensus.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/filterReads/SRR577482.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/filterReads/SRR577483.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/filterReads/GSM1603213.filt.nodup.bam
#samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K4me3_consensus.bam

### RNAseq
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/RNAseq.bam /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/RNAseq_consensus.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/RNAseq.bam.bai /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/RNAseq_consensus.bam.bai

### MNase-seq
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/filterReads/Jkt106_MNase_Seq_20170218_1_val_1Jkt106_MNase_Seq_20170218_2_val_2.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/MNase_single.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/filterReads/Jkt106_MNase_Seq_20170218_1_val_1Jkt106_MNase_Seq_20170218_2_val_2.filt.nodup.bam.bai /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/MNase_single.bam.bai

### DNase-seq
#samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/DNase_consensus.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/filterReads/ENCFF001DPF.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/filterReads/ENCFF001DPG.filt.nodup.bam
#samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/DNase_consensus.bam

### tt-Seq
samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/TTseq_consensus.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam
samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/TTseq_consensus.bam

######################################################################
######################################################################
######################################################################
### Do Histones plus MNase, DNase, RNase, and TT-seq
### make csv file
echo "experiment,file
H3K27ac,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/H3K27ac_consensus.bam
H3K36me3,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/H3K36me3_single.bam
H3K79me3,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/H3K79me3_single.bam
H3K9me3,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/H3K9me3_single.bam
H3K4me1,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/H3K4me1_single.bam
H3K4me3,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/H3K4me3_consensus.bam
RNAse,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/RNAseq_consensus.bam
MNase,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/MNase_single.bam
DNase,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/DNase_consensus.bam
TTseq,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/TTseq_consensus.bam" >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/HIVexpression_HisPlus4.csv

for file in $(ls /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/*0.bed); do
  name=$(basename -s .bed ${file})
  python /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/machine_learning/venkat_rpkm.py --peaks ${file} --experiments /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/HIVexpression_HisPlus4.csv -f ${name}_HisPlus4 --minimum 0
done

######################################################################
######################################################################
######################################################################
### Get chromHMM
grep -v "_" /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/reorder/jurkat_15_segments_reorder.bed | sort -k 1,1 -k2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/chromHMM_histones.bed

for file in $(ls /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/*0.bed); do
  name=$(basename -s .bed ${file})
  bedtools intersect -wo -f 0.5 -a ${file} -b /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/chromHMM_histones.bed  | awk '{OFS="\t"}; {print $1, $2, $3, $4, $8}' >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/${name}_chromHMM.bed
done

######################################################################
######################################################################
######################################################################
### Get lamin; note there are missing values
sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments_hg38.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/laminB_subcompartments.bed

for file in $(ls /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/*0.bed); do
  name=$(basename -s .bed ${file})
  bedtools intersect -wo -f 0.5 -a ${file} -b /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/laminB_subcompartments.bed  | awk '{OFS="\t"}; {print $1, $2, $3, $4, $8}' | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/${name}_lamin.bed
done

######################################################################
######################################################################
######################################################################
# Combine the files and add expression table
awk '{OFS="\t"} {print $1, $2, $3,$4, $5}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.bed | sort -k 1,1 -k2,2n | uniq -u >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/hiv_expression_table_withname.txt
### add expression values to output files; R
module load R
Rscript /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/machine_learning/merge_and_add_expression.R


