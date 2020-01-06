#! /bin/bash

#SBATCH --job-name=LH
#SBATCH --partition=super
#SBATCH --output=LH.%j.out
#SBATCH --error=LH.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load bedtools/2.26.0
module load pybedtools
module load python/2.7.x-anaconda
module load samtools

 
### Get HIV expression in Lamin sub compartment
HIVexp='/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression_sorted.bed'
lamin='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments_hg38.bed'
loop='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/Loop_outercoord.bed'
tad='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/Tad_coord.bed'

#bedtools intersect -wao -a ${HIVexp} -b <(sort -k 1,1 -k 2,2n ${lamin}) | awk '{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/HIVexp_lamin.bed

### Add Loop to file
#bedtools intersect -wao -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/HIVexp_lamin.bed -b ${loop} | awk '{OFS="\t"} NF{NF-=1};1' | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/HIVexp_lamin_loop.bed

### Add Tad to file
#bedtools intersect -wao -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/HIVexp_lamin_loop.bed -b ${tad} | awk '{OFS="\t"} NF{NF-=1};1' | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/HIVexp_lamin_loop_tad.bed

### Get only those in A1, in 1 loop, and in 1 tad
#awk '$10 == "A1" && $11 ~ /chr/ && $14 ~ /chr/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/HIVexp_lamin_loop_tad.bed | sort -u -k 1,3 >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/HIVexp_A1lamin_inloopuniq_intad.bed

### Split the files by compartment type labeling by barcode
#awk '{OFS = "\t"} {print $7, $8, $9, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/HIVexp_A1lamin_inloopuniq_intad.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/Lamin2score.bed
#awk '{OFS = "\t"} {print $11, $12, $13, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/HIVexp_A1lamin_inloopuniq_intad.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/Loop2score.bed
#awk '{OFS = "\t"} {print $14, $15, $16, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/HIVexp_A1lamin_inloopuniq_intad.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/Tad2score.bed

### Calculate RPKM
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
H3K27me3,/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/input_files/H3K27me3_single.bam" >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/score_input.csv


python /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/machine_learning/venkat_rpkm.py --peaks /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/Lamin2score.bed --experiments /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/score_input.csv -f Lamin --minimum 0

python /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/machine_learning/venkat_rpkm.py --peaks /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/Loop2score.bed --experiments /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/score_input.csv -f Loop --minimum 0

python /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/machine_learning/venkat_rpkm.py --peaks /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/Tad2score.bed --experiments /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/lamin_vs_histones/score_input.csv -f Tad --minimum 0
