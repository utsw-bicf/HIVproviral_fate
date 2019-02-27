#! /bin/bash

#SBATCH --job-name=checksum
#SBATCH --partition=super
#SBATCH --output=cs.%j.out
#SBATCH --error=cs.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=NONE

find /archive/shared/DOrso_BICF/ -name "*.fastq.gz" -exec md5sum {} \; >data_checksum_list.txt
