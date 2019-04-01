#! /bin/sh
### Run Hi-C data through Juicer

#SBATCH --job-name=Juicer
#SBATCH --mem 126976
#SBATCH --output=Juicer.%j.out
#SBATCH --error=Juicer.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### Make directory and softlink fastq files
### Read name must be in the form of *_R*.fastq*
#mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/Juicer/fastq
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489136_1.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/Juicer/fastq/GSM3489136_R1.fastq.gz
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489136_2.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/Juicer/fastq/GSM3489136_R2.fastq.gz
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489137_1.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/Juicer/fastq/GSM3489137_R1.fastq.gz
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489137_2.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/Juicer/fastq/GSM3489137_R2.fastq.gz

### Make an enzyme cut file
#module load python/2.7.5
#python /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Chromatin_organization/generate_site_positions.py MboI GRCh38 /project/shared/bicf_workflow_ref/human/GRCh38/genome.fa
#module rm python

module load juicer/1.5.6

juicer.sh -z /project/shared/bicf_workflow_ref/human/GRCh38/genome.fa \
  -p /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes \
  -d /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/Juicer \
  -y /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/Juicer/GRCh38_MboI.txt \
  -q 32GB \
  -l 32GB
  
