#! /bin/bash

########## This script calculates the genome fraction of 
########## the chromHMM states for figure 7A

##### Calculate total count
awk '{sum += $3-$2+1} END {print sum}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed
# 3190548642

##### Loop through each state and calculate percentage
for i in U1 U2 U3 U4 U5 U6 U7 U8 U9 U10 U11 U12 U13 U14 U15; do
  awk -v var="${i}" '{OFS = "\t"} $4 == var {sum += $3-$2+1} END {print var,sum/3190548642*100}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed
done
