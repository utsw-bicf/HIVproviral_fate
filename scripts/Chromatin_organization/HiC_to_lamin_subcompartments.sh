#! /bin/sh
### Run Hi-C data through HOMER

#SBATCH --job-name=lamin
#SBATCH -p 256GB 
#SBATCH --mem 253952
#SBATCH --output=lamin.%j.out
#SBATCH --error=lamin.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load juicebox/1.5.6
module load R
module load bedtools


########################################################################
########################################################################
########################################################################
########## Previously done in HiC_to_lamin_subcompartments.sh
########################################################################
########################################################################
########################################################################

########## Assign Hi-C results to Lamin A/B; and split A/B compartment
########## Combine tag directories
#echo "Combine tag directories"
#makeTagDirectory /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged -d /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1 /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2 -tbp 1

########## Calculate PCA analysis; 
########## Need to do each chromosome separately and flip values on PC
########## In order to keep Lamin A as positive and Lamin B as negative
echo "Calculate PCA"
export TMPDIR=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/tmp
runHiCpca.pl auto \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged \
  -cpu 10

# make chromHMM file
mkdir chromHMM_PC1_input
for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22; do
  grep -w ${i} /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged.50x100kb.PC1.bed | awk '$6 >= 0 {print $0}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_PC1_input/${i}_pos.bed
grep -w ${i} /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged.50x100kb.PC1.bed | awk '$6 < 0 {print $0}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_PC1_input/${i}_neg.bed
done

unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_PC1_input/ \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_PC1_input/test_pos_neg



##### Change PC1 values if Lamin A isn't positive; exclude chrX, chrY, chrM
##### Do this for chromosomes: chr2,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr16,chr17,chr18,chr19,chr20,chr21
awk '{OFS = "\t"} {if ($1 ~/chr2-/ || $1 ~/chr4-/ || $1 ~/chr5-/ || $1 ~/chr6-/ || $1 ~/chr7-/ || $1 ~/chr8-/ || $1 ~/chr9-/ || $1 ~/chr10-/ || $1 ~/chr11-/ || $1 ~/chr12-/ || $1 ~/chr16-/ || $1 ~/chr17-/ || $1 ~/chr18-/ || $1 ~/chr19-/ || $1 ~/chr20-/ || $1 ~/chr21-/) {print $1, $2, $3, $4, $5, $6*-1 } else {print $0}}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged.50x100kb.PC1.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged_fixed.50x100kb.PC1.txt


########## Separate positive (A) and negative (B)
grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged_fixed.50x100kb.PC1.txt | awk '{OFS = "\t"} $6 >= 0 {print $2, $3, $4, "LaminA", $6, $5}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1.bed
grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged_fixed.50x100kb.PC1.txt | awk '{OFS = "\t"} $6 < 0 {print $2, $3, $4, "LaminB", $6, $5}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1.bed
# merge
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/Jurkat_PC1_Lamin.bed

########## Separate LaminA by k-means clustering (2)
Rscript /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Chromatin_organization/kmeans_Lamin.R

########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################



########## Make chromHMM of A1, A2,B vs marks vs GSM
grep "LaminA1" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1_kmeans.bed > /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2/input/jurkat_A1.bed
grep "LaminA2" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1_kmeans.bed > /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2/input/jurkat_A2.bed
cp /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2/input/jurkat_B.bed

grep "A1" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2/input/GSE63525_GM12878_LaminA1.bed
grep "A2" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2/input/GSE63525_GM12878_LaminA2.bed
grep "B1" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2/input/GSE63525_GM12878_LaminB1.bed
grep "B2" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2/input/GSE63525_GM12878_LaminB2.bed
grep "B3" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2/input/GSE63525_GM12878_LaminB3.bed
grep "B4" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2/input/GSE63525_GM12878_LaminB4.bed

# chromHMM
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2/input/ \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2/A1A2B_compare

# jitter




########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################



########## Run Kmeans of 5 on all PC1
Rscript /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Chromatin_organization/kmeans5_allPC1_Lamin.R

########## Make chromHMM of each "K" vs marks vs GSM
grep "A1" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans5/input/GSE63525_GM12878_LaminA1.bed
grep "A2" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans5/input/GSE63525_GM12878_LaminA2.bed
grep "B1" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans5/input/GSE63525_GM12878_LaminB1.bed
grep "B2" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans5/input/GSE63525_GM12878_LaminB2.bed
grep "B3" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans5/input/GSE63525_GM12878_LaminB3.bed
grep "B4" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans5/input/GSE63525_GM12878_LaminB4.bed

awk '{OFS = "\t"} $7 >= 1 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged.50x100kb.PC1_kmeans5.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans5/input/jurkat_K1.bed
awk '{OFS = "\t"} $7 >= 5 {print $1, $2, $3, "K5", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged.50x100kb.PC1_kmeans5.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans5/input/jurkat_K5.bed
awk '{OFS = "\t"} $7 >= 3 {print $1, $2, $3, "K3", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged.50x100kb.PC1_kmeans5.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans5/input/jurkat_K3.bed
awk '{OFS = "\t"} $7 >= 4 {print $1, $2, $3, "K4", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged.50x100kb.PC1_kmeans5.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans5/input/jurkat_K4.bed
awk '{OFS = "\t"} $7 >= 2 {print $1, $2, $3, "K2", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged.50x100kb.PC1_kmeans5.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans5/input/jurkat_K2.bed

unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans5/input/ \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans5/A1A2B1B2B3_compare


########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################


########## Treat A with kmeans 2, and B with kmeans 3
grep "LaminA1" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1_kmeans.bed > /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2_A_kmeans3_B/input/jurkat_A1.bed
grep "LaminA2" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1_kmeans.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2_A_kmeans3_B/input/jurkat_A2.bed
awk '{OFS = "\t"} $7 == 1 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2_A_kmeans3_B/input/LaminB_kmeans_1.bed
awk '{OFS = "\t"} $7 == 2 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2_A_kmeans3_B/input/LaminB_kmeans_2.bed
awk '{OFS = "\t"} $7 == 3 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2_A_kmeans3_B/input/LaminB_kmeans_3.bed

grep "A1" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2_A_kmeans3_B/input/GSE63525_GM12878_LaminA1.bed
grep "A2" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2_A_kmeans3_B/input/GSE63525_GM12878_LaminA2.bed
grep "B1" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2_A_kmeans3_B/input/GSE63525_GM12878_LaminB1.bed
grep "B2" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2_A_kmeans3_B/input/GSE63525_GM12878_LaminB2.bed
grep "B3" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2_A_kmeans3_B/input/GSE63525_GM12878_LaminB3.bed
grep "B4" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2_A_kmeans3_B/input/GSE63525_GM12878_LaminB4.bed


unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2_A_kmeans3_B/input/ \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans2_A_kmeans3_B/A1A2B1B2B3_separatekmeans_compare

########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########## Make supplementary file of A compartment kmeans 2
########## B compartment kmeans 2-5
########## against chromHMM and GSM12878

### Run R script of kmeans
Rscript /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Chromatin_organization/kmeans_Lamin.R


### Separate kmeans
# A
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsA
awk '{OFS = "\t"} $7 == 1 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1_kmeans.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsA/jurkat_subA_kmean1of2.bed
awk '{OFS = "\t"} $7 == 2 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1_kmeans.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsA/jurkat_subA_kmean2of2.bed

# B; kmeans2
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk2
awk '{OFS = "\t"} $7 == 1 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_2.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk2/jurkat_subB_kmean1of2.bed
awk '{OFS = "\t"} $7 == 2 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_2.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk2/jurkat_subB_kmean2of2.bed

# B; kmeans3
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk3
awk '{OFS = "\t"} $7 == 1 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk3/jurkat_subB_kmean1of3.bed
awk '{OFS = "\t"} $7 == 2 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk3/jurkat_subB_kmean2of3.bed
awk '{OFS = "\t"} $7 == 3 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk3/jurkat_subB_kmean3of3.bed

# B; kmeans4
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk4
awk '{OFS = "\t"} $7 == 1 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_4.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk4/jurkat_subB_kmean1of4.bed
awk '{OFS = "\t"} $7 == 2 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_4.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk4/jurkat_subB_kmean2of4.bed
awk '{OFS = "\t"} $7 == 3 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_4.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk4/jurkat_subB_kmean3of4.bed
awk '{OFS = "\t"} $7 == 4 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_4.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk4/jurkat_subB_kmean4of4.bed

# B; kmeans5
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk5
awk '{OFS = "\t"} $7 == 1 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_5.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk5/jurkat_subB_kmean1of5.bed
awk '{OFS = "\t"} $7 == 2 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_5.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk5/jurkat_subB_kmean2of5.bed
awk '{OFS = "\t"} $7 == 3 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_5.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk5/jurkat_subB_kmean3of5.bed
awk '{OFS = "\t"} $7 == 4 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_5.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk5/jurkat_subB_kmean4of5.bed
awk '{OFS = "\t"} $7 == 4 {print $1, $2, $3, "K1", $6, $4}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_5.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk5/jurkat_subB_kmean5of5.bed

# GSM12878
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsGSM12878
grep "A1" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsGSM12878/GM12878_LaminA1.bed
grep "A2" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsGSM12878/GM12878_LaminA2.bed
grep "B1" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsGSM12878/GM12878_LaminB1.bed
grep "B2" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsGSM12878/GM12878_LaminB2.bed
grep "B3" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsGSM12878/GM12878_LaminB3.bed
grep "B4" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsGSM12878/GM12878_LaminB4.bed


### Make chromHMM plots
# A subcomparment, kmeans2
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsA \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/Jurkat_LaminA

# B subcomparment, kmeans2
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk2 \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/Jurkat_LaminB_k2

# B subcomparment, kmeans3
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk3 \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/Jurkat_LaminB_k3

# B subcomparment, kmeans4
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk4 \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/Jurkat_LaminB_k4

# B subcomparment, kmeans5
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsBk5 \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/Jurkat_LaminB_k5

# GSM12878
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/inputsGSM12878 \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_supplementary/GSM12878
