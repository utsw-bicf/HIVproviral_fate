#! /bin/bash

#SBATCH --job-name=TLRM
#SBATCH --partition=super
#SBATCH --output=TLRM.%j.out
#SBATCH --error=TLRM.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load pybedtools
module load python/2.7.x-anaconda
module load samtools
module load bedtools

########## This scores Tads and loops (and not) by rpkm values of MNase seq
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
TTseq,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/TTseq_consensus.bam
H3K27me3,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/H3K27me3_single.bam" >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_TadLoop/MNase_input.csv

echo "Start Tad"
awk '{OFS = "\t"}{print $1, $2, $3}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_TadLoop/Tad.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_TadLoop/Tad_r.bed
python /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/machine_learning/venkat_rpkm_3.py --peaks /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_TadLoop/Tad_r.bed --experiments /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_TadLoop/MNase_input.csv -f tad --minimum 0

echo "Start not Tad"
awk '{OFS = "\t"}{print $1, $2, $3}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/not_tad.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/not_tad_r.bed
python /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/machine_learning/venkat_rpkm_3.py --peaks /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/not_tad_r.bed --experiments /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_TadLoop/MNase_input.csv -f not_tad --minimum 0

echo "Start loop"
awk '{OFS = "\t"}{print $1, $2, $3}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_TadLoop/Loop.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_TadLoop/Loop_r.bed
python /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/machine_learning/venkat_rpkm_3.py --peaks /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_TadLoop/Loop_r.bed --experiments /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_TadLoop/MNase_input.csv -f loop --minimum 0

echo "Start not loop"
awk '{OFS = "\t"}{print $1, $2, $3}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/not_loop.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/not_loop_r.bed
python /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/machine_learning/venkat_rpkm_3.py --peaks /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/not_loop_r.bed --experiments /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_TadLoop/MNase_input.csv -f not_loop --minimum 0
