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

bamCoverage -p 30 --filterRNAstrand reverse -b /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.dedup.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.reverse.bw

### Bedgraph
bamCoverage -p 30 -of bedgraph --filterRNAstrand forward -b /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.forward.bedGraph

bamCoverage -p 30 -of bedgraph --filterRNAstrand reverse -b /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.reverse.bedGraph

bamCoverage -p 30 -of bedgraph --filterRNAstrand forward -b /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.dedup.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.forward.bedGraph

bamCoverage -p 30 -of bedgraph --filterRNAstrand reverse -b /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.dedup.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.reverse.bedGraph


### Bedgraph to bw; note initial bw didn't work for some reason
module load UCSC_userApps/v317
grep "^chr" /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.forward.bedGraph | grep -v "_" | grep -v "chrEBV" >/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.forward_chronly.bedGraph
bedSort /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.forward_chronly.bedGraph /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.forward_sorted.bedGraph
bedGraphToBigWig /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.forward_sorted.bedGraph /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.forward_sorted.bw

grep "^chr" /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.reverse.bedGraph | grep -v "_" | grep -v "chrEBV" >/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.reverse_chronly.bedGraph
bedSort /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.reverse_chronly.bedGraph /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.reverse_sorted.bedGraph
bedGraphToBigWig /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.reverse_sorted.bedGraph /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.reverse_sorted.bw

grep "^chr" /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.forward.bedGraph | grep -v "_" | grep -v "chrEBV" >/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.forward_chronly.bedGraph
bedSort /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.forward_chronly.bedGraph /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.forward_sorted.bedGraph
bedGraphToBigWig /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.forward_sorted.bedGraph /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.forward_sorted.bw

grep "^chr" /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.reverse.bedGraph | grep -v "_" | grep -v "chrEBV" >/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.reverse_chronly.bedGraph
bedSort /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.reverse_chronly.bedGraph /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.reverse_sorted.bedGraph
bedGraphToBigWig /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.reverse_sorted.bedGraph /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.reverse_sorted.bw
