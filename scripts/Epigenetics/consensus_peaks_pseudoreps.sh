#! /bin/bash
### This script takes H3K4 experiments and H3K27ac experiments and calls consensus peaks
### This also tests the difference between using intersectBed and GRanges in R


########## H3K4
### intersectBed

Idir='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/'
H3K4me3="${Idir}SRR061743_H3K4m3_pooled_peaks.narrowPeak ${Idir}SRR577482_H3K4me3_pooled_peaks.narrowPeak ${Idir}SRR577483_H3K4me3_pooled_peaks.narrowPeak ${Idir}GSM1603213_GFP_H3K4me3_pooled_peaks.narrowPeak"

intersectBed -wo -a ${H3K4me3}
  -b %s' % (pool_peaks, true_rep_peaks[0]),
                    awk_command,
                    cut_command,
                    'sort -u']

