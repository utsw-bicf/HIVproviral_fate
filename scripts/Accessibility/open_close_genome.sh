#! /bin/bash

########## This script labeles the genome open/closed
########## This is based on MNase and DNAse

grep -v "log" /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/iNPS/iNPS_output/MNase_significant_all.bed | sort -k 1,1 -k 2,2n | uniq > /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/MNase_peaks.bed

sort -k 1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/consensusPeaks/DNAse-seq_ENCODE_GSM736501.replicated.narrowPeak | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/DNase_peaks.bed

MNpeaks='/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/MNase_peaks.bed'
DNpeaks='/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/DNase_peaks.bed'

### Get regions that MNase only
module load bedtools
bedtools subtract -a ${MNpeaks} -b ${DNpeaks} | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/closed_MNaseonly_peaks.bed

### Get regions that are DNase only
bedtools subtract -a ${DNpeaks} -b ${MNpeaks} | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/open_DNaseonly_peaks.bed

### Get regions that are both
bedtools intersect -a ${DNpeaks} -b ${MNpeaks} | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/both_peaks.bed

### Relabel each and combine
awk '{OFS = "\t"} $1 !~ /_/ {print $1, $2, $3, "closed_MNase"}' /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/closed_MNaseonly_peaks.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/closed_MNaseonly_fixed.bed
awk '{OFS = "\t"} $1 !~ /_/ {print $1, $2, $3, "open_DNase"}' /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/open_DNaseonly_peaks.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/open_DNaseonly_fixed.bed
awk '{OFS = "\t"} $1 !~ /_/ {print $1, $2, $3, "unknown_both"}' /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/both_peaks.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/both_fixed.bed

cat /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/closed_MNaseonly_fixed.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/open_DNaseonly_fixed.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/both_fixed.bed | sort -k 1,1 -V -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/delete.bed


### Find regions that are neither
# fix chrom sizes
sort -k 1,1 -V -k 2,2n /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes | grep -v "_" >/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/chrom.sizes.fixed
bedtools complement -i /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/delete.bed -g /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/chrom.sizes.fixed >/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/neither_peaks.bed

### Add label and combine with previous
awk '{OFS = "\t"} $1 !~ /_/ {print $1, $2, $3, "neither"}' /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/neither_peaks.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/neither_fixed.bed

cat /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/delete.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/neither_fixed.bed | sort -k 1,1 -V -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/label_genome_open_close/jurkat_open_closed.bed
