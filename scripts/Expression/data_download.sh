#! /bin/sh
### Download RNAseq data for d'Orso lab

#SBATCH --job-name=dd
#SBATCH --partition=super
#SBATCH --output=dd.%j.out
#SBATCH --error=dd.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=NONE

module load sra_toolkit/2.8.2-1

### GSM2260195: SRR4000396
### paired-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR400/SRR4000396/SRR4000396.sra
#fastq-dump -I --split-files SRR4000396
#fastq-dump --gzip --skip-technical  --readids --dumpbase --split-files --clip SRR4000396


### GSM2453329: SRR5170980
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR517/SRR5170980/SRR5170980.sra
#fastq-dump --gzip --skip-technical  --readids --dumpbase --split-files --clip SRR5170980
