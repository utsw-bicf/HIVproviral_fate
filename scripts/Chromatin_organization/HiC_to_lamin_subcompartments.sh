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


########## Determine if each chromosome positive value is LaminA or LaminB
##### Copy other Lamin annotations
grep "A1\|A2" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminA.bed
grep "B1\|B2\|B3" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminB.bed

##### Make bed file
grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged.50x100kb.PC1.txt | awk '{OFS = "\t"} {print $2, $3, $4, $5, "NA", $6}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged.50x100kb.PC1.bed

##### Intersect bed file with Lamin to tell direction
bedtools intersect -wao -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged.50x100kb.PC1.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminA.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_PC1_OL_laminA.bed
bedtools intersect -wao -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged.50x100kb.PC1.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminB.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_PC1_OL_laminB.bed

# make txt file
echo -e "chrom\tLaminA\tLaminB" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/Lamin_chromosome_orientation.txt
for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM; do
  LaminA=`grep -w ${i} /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_PC1_OL_laminA.bed | awk '{sum += $16} END {print sum}'`
  LaminB=`grep -w ${i} /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_PC1_OL_laminB.bed | awk '{sum += $16} END {print sum}'`
  echo -e "${i}\t${LaminA}\t${LaminB}" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/Lamin_chromosome_orientation.txt
done

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




########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########## Separate positive (A) and negative (B)
grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged_fixed.50x100kb.PC1.txt | awk '{OFS = "\t"} $6 >= 0 {print $2, $3, $4, "LaminA", $6, $5}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1.bed
grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged_fixed.50x100kb.PC1.txt | awk '{OFS = "\t"} $6 < 0 {print $2, $3, $4, "LaminB", $6, $5}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1.bed
# merge
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/Jurkat_PC1_Lamin.bed

########## Separate LaminA by k-means clustering (2) or
########## Separate LaminB by K-means clustering (3 or 4)
Rscript /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Chromatin_organization/kmeans_Lamin.R

### Combine A1, A2, and B
awk '{OFS = "\t"} {if ($7 == 1) {print $1, $2, $3, "LaminA1", $5, $6} else {print $1, $2, $3, "LaminA2", $5, $6}}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1_kmeans.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1_kmeans.bed
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1_kmeans.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/Jurkat_PC1withKmeans_Lamin.bed

########## Break by K mean and label as A1 or A2 (later), and also B's
awk '{OFS = "\t"} $7 == 2 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1_kmeans.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1_kmeans_2.bed

awk '{OFS = "\t"} $7 == 1 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1_kmeans.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1_kmeans_1.bed

# B, 3
awk '{OFS = "\t"} $7 == 1 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_1of3.bed
awk '{OFS = "\t"} $7 == 2 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_2of3.bed
awk '{OFS = "\t"} $7 == 3 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3of3.bed

# B, 4
awk '{OFS = "\t"} $7 == 1 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_4.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_1of4.bed
awk '{OFS = "\t"} $7 == 2 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_4.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_2of4.bed
awk '{OFS = "\t"} $7 == 3 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_4.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3of4.bed
awk '{OFS = "\t"} $7 == 4 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_4.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_4of4.bed

# B, not chr19
awk '{OFS = "\t"} $7 == 1 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans3_notChr19.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_1_notChr19.bed
awk '{OFS = "\t"} $7 == 2 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans3_notChr19.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_2_notChr19.bed
awk '{OFS = "\t"} $7 == 3 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans3_notChr19.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3_notChr19.bed
########## View in IGV
awk '{OFS = "\t"} $7 == 1 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans4_Chr19.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_1_Chr19.bed
awk '{OFS = "\t"} $7 == 2 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans4_Chr19.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_2_Chr19.bed
awk '{OFS = "\t"} $7 == 3 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans4_Chr19.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3_Chr19.bed
awk '{OFS = "\t"} $7 == 4 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans4_Chr19.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_4_Chr19.bed

# chr19 but with 3 clusters
awk '{OFS = "\t"} $7 == 1 {print $1, $2, $3, $4, $5, $6}' LaminB_negative_PC1_kmeans3_Chr19.txt >chromHMM_kmeans_input/LaminB_19/chr19_1.bed
awk '{OFS = "\t"} $7 == 2 {print $1, $2, $3, $4, $5, $6}' LaminB_negative_PC1_kmeans3_Chr19.txt >chromHMM_kmeans_input/LaminB_19/chr19_2.bed
awk '{OFS = "\t"} $7 == 3 {print $1, $2, $3, $4, $5, $6}' LaminB_negative_PC1_kmeans3_Chr19.txt >chromHMM_kmeans_input/LaminB_19/chr19_3.bed

# Lamin B, Kmeans 3, all chromosomes together
awk '{OFS = "\t"} $7 == 1 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/LaminB_all/LaminB_kmeans_1.bed
awk '{OFS = "\t"} $7 == 2 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/LaminB_all/LaminB_kmeans_2.bed
awk '{OFS = "\t"} $7 == 3 {print $1, $2, $3, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/LaminB_all/LaminB_kmeans_3.bed
# Divide B into chr19 and not, then divide on B1, B2, B3, B4
awk '$4 == "B1" {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminB.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/LaminB_all/GSE63525_LaminB1.bed
awk '$4 == "B2" {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminB.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/LaminB_all/GSE63525_LaminB2.bed
awk '$4 == "B3" {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminB.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/LaminB_all/GSE63525_LaminB3.bed
awk '$4 == "B4" {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminB.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/LaminB_all/GSE63525_LaminB4.bed

########## Compare old vs new in chromHMM
mkdir chromHMM_kmeans_input
### make chromHMM file
# Divide A by A1, A2
awk '$4 == "A1" {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminA.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/GSE63525_LaminA1.bed
awk '$4 == "A2" {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminA.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/GSE63525_LaminA2.bed

# Divide B into chr19 and not, then divide on B1, B2, B3, B4
awk '$4 == "B1" && $1 !~ /chr19/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminB.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/GSE63525_LaminB1_not19.bed
awk '$4 == "B2" && $1 !~ /chr19/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminB.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/GSE63525_LaminB2_not19.bed
awk '$4 == "B3" && $1 !~ /chr19/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminB.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/GSE63525_LaminB3_not19.bed

awk '$4 == "B1" && $1 ~ /chr19/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminB.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/GSE63525_LaminB1_19.bed
awk '$4 == "B2" && $1 ~ /chr19/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminB.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/GSE63525_LaminB2_19.bed
awk '$4 == "B4" && $1 ~ /chr19/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminB4.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/GSE63525_LaminB4_19.bed



unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/LaminA/ \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/test_kmeans_LaminA

unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/LaminB_not19/ \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/test_kmeans_LaminB_not19

unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/LaminB_19/ \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/test_kmeans_LaminB_19


unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/LaminB_all/ \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/chromHMM_kmeans_input/test_kmeans_LaminB_all


##### It appears that Lamin A compartments are ok. And B compartments, not on chr19 are ok.

##### I am skipping B4 and keeping B1, B2, B3
### Combine A1, A2, and B
awk '{OFS = "\t"} {if ($7 == 1) {print $1, $2, $3, "LaminB2", $5, $6} else if ($7 == 2) {print $1, $2, $3, "LaminB3", $5, $6} else {print $1, $2, $3, "LaminB1", $5, $6}}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans_3.txt | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans3.bed
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_positive_PC1_kmeans.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_negative_PC1_kmeans3.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/Jurkat_PC1withKmeans_Lamin_A1A2B1B2B3.bed


########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########## Might not use
########## Find compartments
#echo "Find compartments"
# A region
#findHiCCompartments.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged_fixed.50x100kb.PC1.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged_Acompartments.txt
#sort -k 2,2 -k 3,3n /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged_Acompartments.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged_Acompartments_sorted.txt
## B region
#findHiCCompartments.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged_fixed.50x100kb.PC1.txt -opp >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged_Bcompartments.txt
#sort -k 2,2 -k 3,3n /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged_Bcompartments.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged_Bcompartments_sorted.txt 

########## Make bed file of compartments
#awk '{OFS = "\t"} {print $2, $3, $4, $5, $6, $7}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged_Acompartments.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/Lamin_A_compartments_jurkat.bed
#awk '{OFS = "\t"} {print $2, $3, $4, $5, $6, $7}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged_Bcompartments.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/Lamin_B_compartments_jurkat.bed

########## Test if there is an overlap; Note: we are good to go
#module load bedtools
#bedtools intersect -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/Lamin_A_compartments_jurkat.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/Lamin_B_compartments_jurkat.bed | wc -l






########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########## Find % overlap of A/B compartments between jurkat and "other"
#grep "B4" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminB4.bed

#bedtools intersect -wao -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/Lamin_A_compartments_jurkat.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminA.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_overlap_compare.bed

#overlapA=`awk '{sum += $16} END {print sum}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_overlap_compare.bed`
#LaminA=`awk '{sum += $3-$2} END {print sum}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminA_overlap_compare.bed`
#calc(){ awk "BEGIN { print "$*" }"; }
#calc ${overlapA}/${LaminA}*100

#bedtools intersect -wao -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/Lamin_B_compartments_jurkat.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/GSE63525_GM12878_LaminB.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_overlap_compare.bed
#overlapB=`awk '{sum += $16} END {print sum}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_overlap_compare.bed`
#LaminB=`awk '{sum += $3-$2} END {print sum}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/LaminB_overlap_compare.bed`
#calc(){ awk "BEGIN { print "$*" }"; }
#calc ${overlapB}/${LaminB}*100
