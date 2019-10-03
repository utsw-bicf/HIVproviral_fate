#! /bin/bash

#SBATCH --job-name=jl
#SBATCH --partition=super
#SBATCH --output=jl.%j.out
#SBATCH --error=jl.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END


### This script identifies loops and the CTCF relationship/orientation.
### Find the overlaps with CTCF jaspar

module load bedtools

########################################
### Filter and sort 
grep -v "_" /project/shared/bicf_workflow_ref/dastk_JASPAR_motifs/reformatted_files_human/CTCF.bed | grep "^chr" | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCF.bed

########################################
### Find overlaps with CTCF jaspar overlaps
## Left merge
bedtools intersect -wao -a <(grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D.bed | sort -k 1,1 -k 2,2n) -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCF.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/new_left_loop_jaspar.bed

## Reorder so right side is on left
awk '{OFS = "\t"} {print $4, $5, $6, $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/new_left_loop_jaspar.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input4right_left_loop_jaspar.bed

bedtools intersect -wao -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input4right_left_loop_jaspar.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCF.bed | sort -k 4,4 -k 5,5n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_loop_jaspar.bed

## Has overlaps on both ends
awk '{OFS = "\t"} $10 ~ /chr/ && $20 ~ /chr/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_loop_jaspar.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_both_loop_jaspar.bed
## Has overlaps on 1 end
awk '{OFS = "\t"} ($10 ~ /chr/ && $20 !~ /chr/) || ($20 ~ /chr/ && $10 !~ /chr/) {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_loop_jaspar.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_1end_loop_jaspar.bed
## Has overlaps on neither end
awk '{OFS = "\t"} $10 !~ /chr/ && $20 !~ /chr/ {print $4, $5, $6, $1, $2, $3, $7, $8, $9}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_loop_jaspar.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_noOL_loop_jaspar.bed



########################################
### Find overlaps/Nearest with CTCF jaspar
## Left has overlap
awk '{OFS = "\t"} $10 ~ /chr/ && $20 !~ /chr/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_1end_loop_jaspar.bed | sort -k 1,1 -k 2,2n | cut -f -19 >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCF_1endLeft_loop_jaspar.bed

bedtools closest -D ref -id -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCF_1endLeft_loop_jaspar.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCF.bed | awk '{OFS = "\t"} {print $4, $5, $6, $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, "."}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_LeftOL_rightN_loop_jaspar.bed

## Right has overlap
awk '{OFS = "\t"} $10 !~ /chr/ && $20 ~ /chr/ {print $4, $5, $6, $1, $2, $3, $7, $8, $9, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_1end_loop_jaspar.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCF_1endRight_loop_jaspar.bed

bedtools closest -D ref -iu -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCF_1endRight_loop_jaspar.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCF.bed | awk '{OFS = "\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $20, $21, $22, $23, $24, $25, $26, $27, $28, ".", $10, $11, $12, $13, $14, $15, $16, $17, $18, $19}'  >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_LeftN_rightOL_loop_jaspar.bed


########################################
### Find Nearest (no overlaps on either end) with CTCF jaspar
## Left, downstream of loop
bedtools closest -D ref -iu -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_noOL_loop_jaspar.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCF.bed | awk '{OFS = "\t"} {print $4, $5, $6, $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, "."}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/delete.bed

## Right, upstream of loop
bedtools closest -D ref -id -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/delete.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCF.bed | awk '{OFS = "\t"} {print $4, $5, $6, $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, "."}'  >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_LeftN_rightN_loop_jaspar.bed

rm /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/delete.bed
rm /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/new_left_loop_jaspar.bed
rm /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCF.bed

########################################
## Cat files
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_both_loop_jaspar.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_LeftOL_rightN_loop_jaspar.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_LeftN_rightOL_loop_jaspar.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_LeftN_rightN_loop_jaspar.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_AllPossible_loop_jaspar.bed


## Convergent + on left, - on right
awk '$13 == "+" && $23 == "-" {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_AllPossible_loop_jaspar.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_convergent_loop_jaspar.bed

## Divergent - on left, + on right
awk '$13 == "-" && $23 == "+" {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_AllPossible_loop_jaspar.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_divergent_loop_jaspar.bed

## same on both ends
awk '($13 == "-" && $23 == "-") || ($13 == "+" && $23 == "+") {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_AllPossible_loop_jaspar.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_same_loop_jaspar.bed


### Record possible orientations of loops
## Number of loops: 7571
## Number of loops with both ends O/L CTCF: 1617
## Number of loops with 1 end O/L CTCF: 3500
## Number of loops with 0 end O/L CTCF: 3717
## Convergent: 4329
## Divergent: 885
## Same: 3620


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
#################### Do the same but with peaks
########################################
### Filter and sort 
sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/chipseq_analysis_4CTCF/workflow/output_CTCF/consensusPeaks/CTCF_YOUNG_GSM1689152.replicated.narrowPeak >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCFpeaks.bed

########################################
### Find overlaps with CTCF jaspar overlaps
## Left merge
bedtools intersect -wao -a <(grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D.bed | sort -k 1,1 -k 2,2n) -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCFpeaks.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/new_left_loop_peaks.bed

## Reorder so right side is on left
awk '{OFS = "\t"} {print $4, $5, $6, $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/new_left_loop_peaks.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input4right_left_loop_peaks.bed

bedtools intersect -wao -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input4right_left_loop_peaks.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCFpeaks.bed | sort -k 4,4 -k 5,5n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_loop_peaks.bed

## Has overlaps on both ends
awk '{OFS = "\t"} $10 ~ /chr/ && $21 ~ /chr/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_loop_peaks.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_both_loop_peaks.bed
## Has overlaps on 1 end
awk '{OFS = "\t"} ($10 ~ /chr/ && $21 !~ /chr/) || ($21 ~ /chr/ && $10 !~ /chr/) {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_loop_peaks.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_1end_loop_peaks.bed
## Has overlaps on neither end
awk '{OFS = "\t"} $10 !~ /chr/ && $21 !~ /chr/ {print $4, $5, $6, $1, $2, $3, $7, $8, $9}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_loop_peaks.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_noOL_loop_peaks.bed


########################################
### Find overlaps/Nearest with CTCF peaks
## Left has overlap
awk '{OFS = "\t"} $10 ~ /chr/ && $21 !~ /chr/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_1end_loop_peaks.bed | sort -k 1,1 -k 2,2n | cut -f -20 >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCF_1endLeft_loop_peaks.bed

bedtools closest -D ref -id -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCF_1endLeft_loop_peaks.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCFpeaks.bed | awk '{OFS = "\t"} {print $4, $5, $6, $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_LeftOL_rightN_loop_peaks.bed

## Right has overlap
awk '{OFS = "\t"} $10 !~ /chr/ && $21 ~ /chr/ {print $4, $5, $6, $1, $2, $3, $7, $8, $9, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_1end_loop_peaks.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCF_1endRight_loop_peaks.bed

bedtools closest -D ref -iu -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCF_1endRight_loop_peaks.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCFpeaks.bed | awk '{OFS = "\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}'  >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_LeftN_rightOL_loop_peaks.bed


########################################
### Find Nearest (no overlaps on either end) with CTCF peaks
## Left, downstream of loop
bedtools closest -D ref -iu -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_noOL_loop_peaks.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCFpeaks.bed | awk '{OFS = "\t"} {print $4, $5, $6, $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, "."}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/delete.bed

## Right, upstream of loop
bedtools closest -D ref -id -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/delete.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCFpeaks.bed | awk '{OFS = "\t"} {print $4, $5, $6, $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, "."}'  >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_LeftN_rightN_loop_peaks.bed

rm /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/delete.bed
rm /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/new_left_loop_peaks.bed
rm /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/input_CTCFpeaks.bed

########################################
## Cat files
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_both_loop_peaks.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_LeftOL_rightN_loop_peaks.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_LeftN_rightOL_loop_peaks.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_LeftN_rightN_loop_peaks.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_loop/CTCF_AllPossible_loop_peaks.bed

### Record possible orientations of loops
## Number of loops: 7571
## Number of loops with both ends O/L CTCF: 2343
## Number of loops with 1 end O/L CTCF: 3645
## Number of loops with 0 end O/L CTCF: 3249


