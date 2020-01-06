#!/bin/bash
#SBATCH --job-name=rnaseq_SE
#SBATCH --partition=super
#SBATCH --output=SE.%j.out
#SBATCH --error=SE.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load nextflow/0.24.1-SNAPSHOT

nextflow run workflow/main.nf \
--design '/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/design_SE.txt' \
--fastqs '/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/raw_data/*.gz' \
--genome '/project/shared/bicf_workflow_ref/GRCh38/' \
--pairs 'se' \
--output 'workflow/output_SE' \
--dea 'skip'

