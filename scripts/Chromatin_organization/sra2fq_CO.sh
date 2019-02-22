#! /bin/sh
### Download HiC, DamID, CTCF data for d'Orso lab
### Chipseq

#SBATCH --job-name=sra2fq
#SBATCH --partition=super
#SBATCH --output=sra2fq.%j.out
#SBATCH --error=sra2fq.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load sra_toolkit/2.8.2-1

### Hi-C, PE
#fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8244644
#fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8244645
#fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8244646
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8388751
rm /home2/s185797/ncbi/public/sra/SRR8388751.sra
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8388752
rm /home2/s185797/ncbi/public/sra/SRR8388752.sra
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8388753
rm /home2/s185797/ncbi/public/sra/SRR8388753.sra 
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8388754
rm /home2/s185797/ncbi/public/sra/SRR8388754.sra
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8388755
rm /home2/s185797/ncbi/public/sra/SRR8388755.sra
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8388756
rm /home2/s185797/ncbi/public/sra/SRR8388756.sra 
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8388757
rm /home2/s185797/ncbi/public/sra/SRR8388757.sra 

### DamID, SE
fastq-dump --gzip --dumpbase --origfmt SRR5261759
rm /home2/s185797/ncbi/public/sra/SRR5261759.sra 
fastq-dump --gzip --dumpbase --origfmt SRR5261760 
rm /home2/s185797/ncbi/public/sra/SRR5261760.sra 
fastq-dump --gzip --dumpbase --origfmt SRR5261761
rm /home2/s185797/ncbi/public/sra/SRR5261761.sra 
fastq-dump --gzip --dumpbase --origfmt SRR5261762 
rm /home2/s185797/ncbi/public/sra/SRR5261762.sra 

### CTCF, SE
# Input
fastq-dump --gzip --dumpbase --origfmt SRR2029842
rm /home2/s185797/ncbi/public/sra/SRR2029842.sra 
# Experiment
fastq-dump --gzip --dumpbase --origfmt SRR2029843
rm /home2/s185797/ncbi/public/sra/SRR2029843.sra

# Experiment, no input given
fastq-dump --gzip --dumpbase --origfmt SRR014986
rm /home2/s185797/ncbi/public/sra/SRR014986.sra
