#! /bin/sh

#SBATCH --job-name=integration
#SBATCH --partition=super
#SBATCH --output=integration.%j.out
#SBATCH --error=integration.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load nextflow/0.24.1-SNAPSHOT
module load singularity/2.6.1
module load python/2.7.5
module load java/oracle/jdk1.8.0_65

nextflow run /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/map.nf \
  --index /project/shared/bicf_workflow_ref/human/GRCh38/genome.fa \
  -with-singularity /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/bioinformatics.simg
