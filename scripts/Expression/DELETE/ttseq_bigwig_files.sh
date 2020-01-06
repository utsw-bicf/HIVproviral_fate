#! /bin/bash
#SBATCH --job-name=ttBW
#SBATCH -p 256GB
#SBATCH --mem 252928
#SBATCH --output=ttBW.%j.out
#SBATCH --error=ttBW.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### This script merges bam files, then converts the bam to bigwig files

module load samtools
module load deeptools/2.5.0.1

### Merge bam files for all, fwd,rev
# fwd
samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_fwd.bam \
 /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.fwd.bam \
 /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.fwd.bam \

samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_fwd.bam

bamCoverage --numberOfProcessors=max/2 \
  -b /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_fwd.bam \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_fwd.bw

# rev
samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_rev.bam \
 /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.rev.bam \
 /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.rev.bam \

samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_rev.bam

bamCoverage --numberOfProcessors=max/2 \
  -b /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_rev.bam \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_rev.bw

# all
bamCoverage --numberOfProcessors=max/2 \
  -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Ttseq_CRAMER_merged.bam \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/TTseq_all.bw


