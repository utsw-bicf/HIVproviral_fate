#! /bin/bash
#### This script makes the 3D folding of the genomes
########## Do Not run this, do each step separately

########## 1) Run Hi-C pro
##### Move fastq files to rep1 and rep2
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489136_1.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/data/rep1/
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489136_2.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/data/rep1/
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489137_1.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/data/rep2/
ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489137_2.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/data/rep2/

##### Make config file
### On 256 web viz

module load singularity
singularity exec --contain --cleanenv -B /project /project/apps/singularity-images/hicpro/hicpro-2.11.3.simg HiC-Pro -i /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/data -o /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419 -c /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/config-hicpro.txt

########## Sum the 2 matrix files in R
##### Do this for the raw files, for each of the 
Rscript /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Chromatin_organization/combine_matrix_genome_folding.R

for x in 10000,20000,40000,150000,500000,1000000 {
  cp /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/rep1/raw/${x}/rep1_${x}_abs.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/${x}/combine_${x}_abs.bed
}

########## Run Pastis; just an example
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/pastis/10000

cd /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/pastis/10000

module load python/2.7.x-anaconda
source activate /work/BICF/s185797/programs/pastis/pastis/conda_environ

echo "[all]
output_name: HiC_combine_10000
verbose: 1
max_iter: 100
counts: /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/10000/combine_10000.matrix
lengths: /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/10000/combine_10000_abs.bed
normalize: True" > /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/pastis/10000/config.ini

pastis-mds .
#pastis-pm1 .
#pastis-nmds .
#pastis-pm2 .

source deactivate



