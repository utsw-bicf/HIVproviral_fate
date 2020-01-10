#! /bin/bash
### Download the singularity container
### Following the instructions in https://www.researchgate.net/publication/324548281_Using_Barcoded_HIV_Ensembles_B-HIVE_for_Single_Provirus_Transcriptomics
### See AS6294129323663361527075145563_content_1.pdf

module load singularity/2.6.1
singularity pull docker://ezorita/bioinformatics /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/

git clone https://github.com/gui11aume/BHIVE_for_single_provirus_transcriptomics.git /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/

### Edit the following files
# map.nf
# map.cfg
# expr.nf
# expr.cfg

