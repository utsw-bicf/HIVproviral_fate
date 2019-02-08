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
#fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR4000396

### GSM2453329: SRR5170980
### single-end
#fastq-dump --gzip --dumpbase --origfmt SRR5170980

### TTseq
### paired-end
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR4000388
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR4000389

### 4sUseq
### single-end
fastq-dump --gzip --dumpbase --origfmt SRR2130286
fastq-dump --gzip --dumpbase --origfmt SRR2130287
fastq-dump --gzip --dumpbase --origfmt SRR2130288
