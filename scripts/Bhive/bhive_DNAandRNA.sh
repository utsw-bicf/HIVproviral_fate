#! /bin/sh

#SBATCH --job-name=DandR
#SBATCH --partition=super
#SBATCH --mem 123904
#SBATCH --output=DandR.%j.out
#SBATCH --error=DandR.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load nextflow/0.24.1-SNAPSHOT
module load singularity/2.6.1
module load python/2.7.5
module load java/oracle/jdk1.8.0_65
module load R/3.3.2-gccmkl

nextflow run /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/expr.nf \
  --integs /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_integrations.txt \
  -with-singularity /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/bioinformatics.simg
