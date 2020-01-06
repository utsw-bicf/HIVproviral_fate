#! /bin/bash

######### This script labels each loop/tad into a chromHMM slot
######### Then plot the HIV expression into these groups

module load bedtools

## Fix loop
grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D.bed | awk '{OFS = "\t"}{print $1, $2, $6, $7, $8, $9}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/loop.bed
cp /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/loop.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/chromHMM_overlapEnrichment/input_db/

## Fix tad
grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.tad.2D.bed | awk '{OFS = "\t"}{print $1, $2, $6, $7, $8, $9}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/tad.bed
cp /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/tad.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/chromHMM_overlapEnrichment/input_db/

## Fix subcompartments
for i in A1 A2 B1 B2 B3 B4; do
  grep -w ${i} /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/subcompartments/GSE63525_GM12878_subcompartments_hg38.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/chromHMM_overlapEnrichment/input_db/lamin_${i}.bed
done

########## Make chromHMM overlap of loops and tads, with lamin
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/chromHMM_overlapEnrichment/input_db \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/chromHMM_overlapEnrichment/overlap_dbs_LoopTad


########## Loop
##### Classify loops into chromHMM
## Fix chromHMM by each U
for i in U1 U2 U3 U4 U5 U6 U7 U8 U9 U10 U11 U12 U13 U14 U15; do
  grep -w ${i} /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/chromHMM_${i}.bed
done

for i in U1 U2 U3 U4 U5 U6 U7 U8 U9 U10 U11 U12 U13 U14 U15; do
  bedtools intersect -wao -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/loop.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/chromHMM_${i}.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/loop_annotated_${i}.bed
done

## Sum up by chromHMM
for i in U1 U2 U3 U4 U5 U6 U7 U8 U9 U10 U11 U12 U13 U14 U15; do
awk -F '\t' '{OFS = "\t"}{arr[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$10]+=$11;} END {for (i in arr) print i, arr[i]}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/loop_annotated_${i}.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/loop_annotated_${i}_collapsed.bed
done

## Merge files together
I='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM'
paste ${I}/loop_annotated_U1_collapsed.bed ${I}/loop_annotated_U2_collapsed.bed ${I}/loop_annotated_U3_collapsed.bed ${I}/loop_annotated_U4_collapsed.bed ${I}/loop_annotated_U5_collapsed.bed ${I}/loop_annotated_U6_collapsed.bed ${I}/loop_annotated_U7_collapsed.bed ${I}/loop_annotated_U8_collapsed.bed ${I}/loop_annotated_U9_collapsed.bed ${I}/loop_annotated_U10_collapsed.bed ${I}/loop_annotated_U11_collapsed.bed ${I}/loop_annotated_U12_collapsed.bed ${I}/loop_annotated_U13_collapsed.bed ${I}/loop_annotated_U14_collapsed.bed ${I}/loop_annotated_U15_collapsed.bed | awk '{OFS ="\t"}{print $1, $2, $3, $4, $5, $6, $8, $16, $24, $32, $40, $48, $56, $64, $72, $80, $88, $96, $104, $112, $120}' >${I}/loop_annotated_top.bed

header='chrom\tstart\tend\tcolor\tinfo\tinfo\tU1\tU2\tU3\tU4\tU5\tU6\tU7\tU8\tU9\tU10\tU11\tU12\tU13\tU14\tU15'
sed -i '1i'${header} ${I}/loop_annotated_top.bed


## Get the chromHMM that's the max
## Do in excel

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
########## Tads
##### Classify tads into chromHMM
for i in U1 U2 U3 U4 U5 U6 U7 U8 U9 U10 U11 U12 U13 U14 U15; do
  bedtools intersect -wao -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/tad.bed -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/chromHMM_${i}.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/tad_annotated_${i}.bed
done

## Sum up by chromHMM
for i in U1 U2 U3 U4 U5 U6 U7 U8 U9 U10 U11 U12 U13 U14 U15; do
awk -F '\t' '{OFS = "\t"}{arr[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$10]+=$11;} END {for (i in arr) print i, arr[i]}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/tad_annotated_${i}.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM/tad_annotated_${i}_collapsed.bed
done

## Merge files together
I='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_tad_chromHMM'
paste ${I}/tad_annotated_U1_collapsed.bed ${I}/tad_annotated_U2_collapsed.bed ${I}/tad_annotated_U3_collapsed.bed ${I}/tad_annotated_U4_collapsed.bed ${I}/tad_annotated_U5_collapsed.bed ${I}/tad_annotated_U6_collapsed.bed ${I}/tad_annotated_U7_collapsed.bed ${I}/tad_annotated_U8_collapsed.bed ${I}/tad_annotated_U9_collapsed.bed ${I}/tad_annotated_U10_collapsed.bed ${I}/tad_annotated_U11_collapsed.bed ${I}/tad_annotated_U12_collapsed.bed ${I}/tad_annotated_U13_collapsed.bed ${I}/tad_annotated_U14_collapsed.bed ${I}/tad_annotated_U15_collapsed.bed | awk '{OFS ="\t"}{print $1, $2, $3, $4, $5, $6, $8, $16, $24, $32, $40, $48, $56, $64, $72, $80, $88, $96, $104, $112, $120}' >${I}/tad_annotated_top.bed

header='chrom\tstart\tend\tcolor\tinfo\tinfo\tU1\tU2\tU3\tU4\tU5\tU6\tU7\tU8\tU9\tU10\tU11\tU12\tU13\tU14\tU15'
sed -i '1i'${header} ${I}/tad_annotated_top.bed


## Get the chromHMM that's the max
## Do in excel


##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
### Plot in R

