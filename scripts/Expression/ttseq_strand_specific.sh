#!/bin/bash

#SBATCH --job-name=strand
#SBATCH --partition=super
#SBATCH --output=strand.%j.out
#SBATCH --error=strand.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load deeptools/2.5.0.1

#bamCoverage -p 30 --filterRNAstrand forward -b /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.forward.bw

#bamCoverage -p 30 --filterRNAstrand reverse -b /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.reverse.bw

#bamCoverage -p 30 --filterRNAstrand forward -b /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.dedup.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.forward.bw

#bamCoverage -p 30 --filterRNAstrand reverse -b /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.dedup.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.reverse.bw

### Bedgraph
#bamCoverage -p 30 -of bedgraph --filterRNAstrand forward -b /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.forward.bedGraph

#bamCoverage -p 30 -of bedgraph --filterRNAstrand reverse -b /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.reverse.bedGraph

#bamCoverage -p 30 -of bedgraph --filterRNAstrand forward -b /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.dedup.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.forward.bedGraph

#bamCoverage -p 30 -of bedgraph --filterRNAstrand reverse -b /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.dedup.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.reverse.bedGraph


### Bedgraph to bw; note initial bw didn't work for some reason
#module load UCSC_userApps/v317
#grep "^chr" /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.forward.bedGraph | grep -v "_" | grep -v "chrEBV" >/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.forward_chronly.bedGraph
#bedSort /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.forward_chronly.bedGraph /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.forward_sorted.bedGraph
#bedGraphToBigWig /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.forward_sorted.bedGraph /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.forward_sorted.bw

#grep "^chr" /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.reverse.bedGraph | grep -v "_" | grep -v "chrEBV" >/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.reverse_chronly.bedGraph
#bedSort /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.reverse_chronly.bedGraph /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.reverse_sorted.bedGraph
#bedGraphToBigWig /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.reverse_sorted.bedGraph /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.reverse_sorted.bw

#grep "^chr" /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.forward.bedGraph | grep -v "_" | grep -v "chrEBV" >/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.forward_chronly.bedGraph
#bedSort /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.forward_chronly.bedGraph /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.forward_sorted.bedGraph
#bedGraphToBigWig /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.forward_sorted.bedGraph /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.forward_sorted.bw

#grep "^chr" /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.reverse.bedGraph | grep -v "_" | grep -v "chrEBV" >/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.reverse_chronly.bedGraph
#bedSort /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.reverse_chronly.bedGraph /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.reverse_sorted.bedGraph
#bedGraphToBigWig /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.reverse_sorted.bedGraph /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.reverse_sorted.bw


#################################################################
#################################################################
#################################################################
### Split the bam files forward and reverse
module load samtools

Floc='/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT'
ttfiles='Ttseq_CRAMER_GSM2260188.dedup.bam Ttseq_CRAMER_GSM2260187.dedup.bam'

for file in ${ttfiles}; do
  name=$(basename ${file} .dedup.bam)
  ### forward
  # include reads that are 2nd in a pair (128);
  # exclude reads that are mapped to the reverse strand (16)
  samtools view -b -f 128 -F 16 ${Floc}/${file} >${Floc}/${name}.fwd1.bam

  # exclude reads that are mapped to the reverse strand (16) and
  # first in a pair (64): 64 + 16 = 80
#  samtools view -b -f 80 ${Floc}/${file} >${Floc}/${name}.fwd2.bam

  # combine the temporary files
  samtools merge -f ${Floc}/${name}.fwd.bam ${Floc}/${name}.fwd1.bam ${Floc}/${name}.fwd2.bam

  # index the filtered BAM file
  samtools index ${Floc}/${name}.fwd.bam

  # run bamCoverage
  bamCoverage -b ${Floc}/${name}.fwd.bam -o ${Floc}/${name}.fwd.bigWig

  # remove the temporary files
  #rm ${name}.fwd*.bam

  ### Reverse
  # include reads that map to the reverse strand (128)
  # and are second in a pair (16): 128 + 16 = 144
#  samtools view -b -f 144 ${Floc}/${file} >${Floc}/${name}.rev1.bam

  # include reads that are first in a pair (64), but
  # exclude those ones that map to the reverse strand (16)
#  samtools view -b -f 64 -F 16 ${Floc}/${file} >${Floc}/${name}.rev2.bam

  # merge the temporary files
#  samtools merge -f ${Floc}/${name}.rev.bam ${Floc}/${name}.rev1.bam ${Floc}/${name}.rev2.bam

  # index the merged, filtered BAM file
#  samtools index ${Floc}/${name}.rev.bam

  # run bamCoverage
  bamCoverage -b ${Floc}/${name}.rev.bam -o ${Floc}/${name}.rev.bw

  # remove temporary files
  #rm a.rev*.bam
done



