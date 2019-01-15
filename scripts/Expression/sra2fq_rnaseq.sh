#! /bin/sh
### Download RNAseq data for d'Orso lab
### Chipseq

#SBATCH --job-name=sra2fq
#SBATCH --partition=super
#SBATCH --output=sra2fq.%j.out
#SBATCH --error=sra2fq.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load sra_toolkit/2.8.2-1

### GSM2260195: SRR4000396
### paired-end
#fastq-dump --gzip --skip-technical  --readids --dumpbase --split-files --clip SRR4000396
#fastq-dump --gzip --dumpbase --origfmt ${sra}
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR4000396

### GSM2453329: SRR5170980
### single-end
#wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR517/SRR5170980/SRR5170980.sra
#fastq-dump --gzip --skip-technical  --readids --dumpbase --split-files --clip SRR5170980
fastq-dump --gzip --dumpbase --origfmt SRR5170980
