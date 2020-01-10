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

### Make HIV bed file
### Need to remove duplicates which are 1)chr3	48247567	48251569, 2)chr4	79730674	79734676
#awk '{OFS="\t"} {print $1, $2-2001, $3+2000}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.bed | sort -k 1,1 -k2,2n | uniq -u >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/hiv_expression_ML.bed

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

### make csv file
echo "experiment,file
H3K27ac,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K27ac_consensus.bam
H3K36me3,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K36me3_single.bam
H3K79me3,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K79me3_single.bam
H3K9me3,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K9me3_single.bam
H3K4me1,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K4me1_single.bam
H3K4me3,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K4me3_consensus.bam" >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/HIVexpression_histones.csv


### RPKM Filtering
### Get RPKM
python /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/machine_learning/venkat_rpkm.py --peaks /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/hiv_expression_ML.bed --experiments /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/HIVexpression_histones.csv -f HIVexpression_histones --minimum 0

######################################################################
######################################################################
######################################################################
### Do Histones plus MNase, DNase, RNase

### RNAseq
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/RNAseq.bam /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/RNAseq_consensus.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/RNAseq.bam.bai /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/RNAseq_consensus.bam.bai

### MNase-seq
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/filterReads/Jkt106_MNase_Seq_20170218_1_val_1Jkt106_MNase_Seq_20170218_2_val_2.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/MNase_single.bam
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/filterReads/Jkt106_MNase_Seq_20170218_1_val_1Jkt106_MNase_Seq_20170218_2_val_2.filt.nodup.bam.bai /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/MNase_single.bam.bai

### DNase-seq
#samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/DNase_consensus.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/filterReads/ENCFF001DPF.filt.nodup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/filterReads/ENCFF001DPG.filt.nodup.bam
#samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/DNase_consensus.bam

### make csv file
echo "experiment,file
H3K27ac,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K27ac_consensus.bam
H3K36me3,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K36me3_single.bam
H3K79me3,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K79me3_single.bam
H3K9me3,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K9me3_single.bam
H3K4me1,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K4me1_single.bam
H3K4me3,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/H3K4me3_consensus.bam
RNAse,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/RNAseq_consensus.bam
MNase,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/MNase_single.bam
DNase,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/DNase_consensus.bam" >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/HIVexpression_HisPlus3.csv

python /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/machine_learning/venkat_rpkm.py --peaks /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/hiv_expression_ML.bed --experiments /project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/HIVexpression_HisPlus3.csv -f HIVexpression_HisPlus3 --minimum 0



######################################################################
######################################################################
######################################################################
# Make expression table
awk '{OFS="\t"} {print $1, $2-2001, $3+2000,$5}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.bed | sort -k 1,1 -k2,2n | uniq -u >/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/hiv_expression_table.txt
### add expression values to output files; R
module load R
Rscript /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/machine_learning/add_expression.R
