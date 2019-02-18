#! /bin/sh
### Step 1: barcodes

#SBATCH --job-name=rtg
#SBATCH --partition=super
#SBATCH --output=rtg.%j.out
#SBATCH --error=rtg.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load nextflow/0.24.1-SNAPSHOT
module load singularity/2.6.1
module load python/2.7.5
module load java/oracle/jdk1.8.0_65
module load bwa/intel/0.7.15
source ~/.bash_rc

nextflow run /work/BICF/s185797/programs/BHIVE_for_single_provirus_transcriptomics-master/map.nf \
  --index /project/shared/bicf_workflow_ref/human/GRCh38/genome.fa \
#  -with-docker ezorita/bioinformatics \
  -resume
