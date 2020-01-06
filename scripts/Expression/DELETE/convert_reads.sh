#! /bin/bash
### This file converts bam file to tag align file
### Input: dedup.bam

module load python/3.6.1-2-anaconda
module load samtools/1.4.1
module load bedtools/2.26.0

### for macs
module load macs/2.1.0-20151222
module load phantompeakqualtools/1.2
module load UCSC_userApps/v317
module load bedtools/2.26.0

module load deeptools/2.5.0.1
module load samtools

### Idea behind script from:
#python3 /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/scripts/convert_reads.py -b /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/work/c0/c46208c6356facd56dcae1ed28b337/SRR2043614.filt.nodup.bam

#bam_files=`ls /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/*.dedup.bam`
bam_files=`ls /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam`

### Use bedtools to convert to tagAlign
for bam in ${bam_files}; do
  #echo ${bam}
  #prefix=$(basename ${bam} .dedup.bam)
  #echo ${prefix}
  #bamToBed -i $bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}'| gzip >/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/convert_reads/${prefix}.tagAlign.gz
###Works
  #python3 /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/scripts/convert_reads.py -b ${bam} -p
  #python3 /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/scripts/xcor.py -t Ttseq_CRAMER_GSM2260187.dedup.tagAlign.gz -p
  #python3 /home2/s185797/Desktop/Holly_git/pipeline_branch_issues/ATACseq/atacseq_analysis/workflow/scripts/call_peaks_macs.py -t Ttseq_CRAMER_GSM2260187.dedup.tagAlign.gz -x Ttseq_CRAMER_GSM2260187.dedup.tagAlign.15.tagAlign.gz -s Ttseq_CRAMER_GSM2260187 -g hs -z /project/shared/bicf_workflow_ref/GRCh38/genomefile.txt -p -a
###

### Try binning and counting coverage
  #samtools index $bam
  #bamCoverage -b $bam --binSize 200 -o coverage.bw
  bamCoverage -b $bam --extendReads --binSize 200 -o coverage2.bw
  bamCoverage -b ../Ttseq_CRAMER_GSM2260187.dedup.bam --binSize 200 --outFileFormat bedgraph -o coverage.bedgraph
  sort -k1,1 -k2,2n coverage.bedgraph >coverage.bed
  bedtools merge -d 200 -i coverage.bed >coverage_out.bed ### doesn't work





done
