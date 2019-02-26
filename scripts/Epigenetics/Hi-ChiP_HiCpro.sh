#! /bin/sh
### This script processes the Hi-ChiP data with HiCpro

#SBATCH --job-name=HiChiP
#SBATCH --partition=super
#SBATCH --output=HiChiP.%j.out
#SBATCH --error=HiChiP.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load hicpro/2.11.1
module load bowtie2/2.3.2
module load samtools/1.6
module load R/3.3.2-gccmkl
module load python/2.7.x-anaconda

### Need 3 files 
### 1. Bed file of restriction fragments after digestion
python /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Epigenetics/digest_genome.py -r MboI -o /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/GRCh39_MbolI.bed /project/shared/bicf_workflow_ref/human/GRCh38/genome.fa

### 2. Chromosome sizes file; but remove scaffolds
# /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes

### 3. Bowtie index
# /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/bowtie2_index
# OR
# /project/shared/bicf_workflow_ref/human/GRCh38/bowtie2_idx

### Organize the files in proper structure

### Prep config-hicpro.txt file in local folder

###Then run
singularity exec /project/apps/singularity-images/hicpro/hicpro-2.11.1.simg HiC-Pro -i /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/data/ -o /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/output -c /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Epigenetics/config-hicpro.txt
