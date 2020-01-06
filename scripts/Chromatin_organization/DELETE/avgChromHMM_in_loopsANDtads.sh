#! /bin/bash
module load bedtools

########## This script will calculate the aveage ChromHMM score for 3 files.
########## 1) HIV insertions, 2) No HIV insertion but has loop/tad, 3) rest of genome
########## Plot violin plot of score

### loop/tad: chr1	start1	end1	chr2	start2	end2	info	info	info
### chromHMM: (jurkat_15_dense.bed, histone, learn_states_15) chr10	0	73800	10	0	.	0	73800	204,153,255

### Fix loop file, to center on each end, and only loops on same chromosome
grep -v "_" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D.bed | awk '{OFS="\t"}; $1 == $4 {print $1, ($3+$2)/2, ($5+$6)/2, $9}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D_filtered.bed

grep -v "_" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.tad.2D.bed | awk '{OFS="\t"}; $1 == $4 {print $1, $2, $3, $8}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.tad.2D_filtered.bed

### Fix chromHMM file
grep -v "_" /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/jurkat_15_dense.bed | awk '{OFS="\t"}; {print $1, $2, $3, $4}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/chromHMM_15_histone.bed

### Create a perl program that calculates the score
perl /Users/holly/Desktop/holly_git/collaborations/issue141_DOrsoIvan/scripts/Chromatin_organization/avgChromHMM_in_loopsANDtads.pl \
  /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/chromHMM_15_histone.bed \
  /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D_filtered.bed \
  >/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D_filtered_scored.bed

perl /Users/holly/Desktop/holly_git/collaborations/issue141_DOrsoIvan/scripts/Chromatin_organization/avgChromHMM_in_loopsANDtads.pl \
  /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/chromHMM_15_histone.bed \
  /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.tad.2D_filtered.bed \
  >/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.tad.2D_filtered_scored.bed

### Find loops/tads that contain or don't contain the HIV expression (insertion)
bedtools intersect -a lib.loop.2D_filtered_scored.bed \
  -b HIV_expression.bed \
  -wo \
  >HIVexpression_scored_loop.bed

bedtools intersect -a lib.tad.2D_filtered_scored.bed \
  -b HIV_expression.bed \
  -wo \
  >HIVexpression_scored_tad.bed

bedtools intersect -a lib.loop.2D_filtered_scored.bed \
  -b HIV_expression.bed \
  -v \
  >noHIV_in_scored_loop.

bedtools intersect -a lib.tad.2D_filtered_scored.bed \
  -b HIV_expression.bed \
  -v \
  >noHIV_in_scored_tad.bed

### Plot average score in each group in R

