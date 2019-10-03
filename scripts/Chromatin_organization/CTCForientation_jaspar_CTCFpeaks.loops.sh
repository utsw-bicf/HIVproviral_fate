#! /bin/bash

########## This program tries to identify the CTCF, loop orientation
########## First use CTCF peaks against possible CTCF jaspar orientations
########## Then annotate the loops

module load bedtools

### Compare jaspar and peaks
grep -v "_" /project/shared/bicf_workflow_ref/dastk_JASPAR_motifs/reformatted_files_human/CTCF.bed | grep "^chr" | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/input_CTCF_jaspar.bed
sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/chipseq_analysis_4CTCF/workflow/output_CTCF/consensusPeaks/CTCF_YOUNG_GSM1689152.replicated.narrowPeak >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/input_CTCFpeaks.bed

bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/input_CTCF_jaspar.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/input_CTCFpeaks.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/input_jaspar_CTCFpeaks.bed


### Find overlaps with CTCF jaspar/peaks overlaps
## Left merge; and reorder for right merge
bedtools intersect -wao -a <(grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D.bed | sort -k 1,1 -k 2,2n) -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/input_jaspar_CTCFpeaks.bed | awk '{OFS = "\t"} {print $4, $5, $6, $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/delete.bed

bedtools intersect -wao -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/delete.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/input_jaspar_CTCFpeaks.bed | awk '{OFS = "\t"} {print $4, $5, $6, $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks.txt


### Find percentage of loops that have a CTCF binding spot (1, or 2+) on both ends
# 635 are uniq loops for both
awk '$10 ~ /chr/ && $20 ~ /chr/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends.txt
awk '$10 !~ /chr/ && $20 !~ /chr/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_neither_ends.txt
awk '($10 ~ /chr/ && $20 !~ /chr/) || ($10 !~ /chr/ && $20 ~ /chr/) {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_one_ends.txt


### Of those with both ends how many only have 1 possibility at both ends
# 468 uniq at each end
awk '{print $8}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends.txt | sort -k 1,1n | uniq -u >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/uniq_both_ends.txt
grep -vFf /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/uniq_both_ends.txt /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends_unique.txt

## More than 1
# 404 of which there are 164 loops (duplicated on average 2.5 X)
grep -vFf /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/uniq_both_ends.txt /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends_duplicates.txt


## Stats; CTCF on both ends and are not duplicated; Convergent, Divergent, same
## Convergent + on left, - on right
awk '$13 == "+" && $23 == "-" {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends_unique.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends_unique_Convergent.txt
## Divergent - on left, + on right
awk '$13 == "-" && $23 == "+" {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends_unique.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends_unique_Divergent.txt
## same on both ends
awk '($13 == "-" && $23 == "-") || ($13 == "+" && $23 == "+") {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends_unique.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends_unique_Same.txt

## Convergent + on left, - on right
awk '$13 == "+" && $23 == "-" {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends_duplicates.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends_duplicates_Convergent.txt
## Divergent - on left, + on right
awk '$13 == "-" && $23 == "+" {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends_duplicates.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends_duplicates_Divergent.txt
## same on both ends
awk '($13 == "-" && $23 == "-") || ($13 == "+" && $23 == "+") {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends_duplicates.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/CTCF_jaspar_peaks_loop/CTCF_loop_jaspar_peaks_both_ends_duplicates_Same.txt

# Loops with CTCF both, uniq: 468
# Loops with CTCF both, dup: 404
# Loops with CTCF 1end: 2913
# Loops with CTCF none: 4391
# CTCF both, uniq, Convergent: 405
# CTCF both, uniq, Divergent: 2
# CTCF both, uniq, Same: 61
# CTCF both, dup, Convergent: 296
# CTCF both, dup, Divergent: 7
# CTCF both, dup, Same: 101


