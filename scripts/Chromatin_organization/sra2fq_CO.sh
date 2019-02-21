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
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8244644
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8244645
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8244646
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8388751
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8388752
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8388753 
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8388754
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8388755
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8388756 
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR8388757 

### DamID, SE
fastq-dump --gzip --dumpbase --origfmt SRR5261759
fastq-dump --gzip --dumpbase --origfmt SRR5261760 
fastq-dump --gzip --dumpbase --origfmt SRR5261761
fastq-dump --gzip --dumpbase --origfmt SRR5261762 

### CTCF, SE
# Input
fastq-dump --gzip --dumpbase --origfmt SRR2029842
# Experiment
fastq-dump --gzip --dumpbase --origfmt SRR2029843

# Experiment, no input given
fastq-dump --gzip --dumpbase --origfmt SRR014986
