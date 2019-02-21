#! /bin/sh
### Download RNAseq data for d'Orso lab
### Chipseq

#SBATCH --job-name=sra2fq
#SBATCH --partition=super
#SBATCH --output=sra2fq.%j.out
#SBATCH --error=sra2fq.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### Run inside directory: /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/
module load sra_toolkit/2.8.2-1

### Bhive
### PE Integration map, rep 1 and 2
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR3614675
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR3614676

## SE DNA barcodes, rep 1 and 2
fastq-dump --gzip --dumpbase --origfmt SRR3614677
fastq-dump --gzip --dumpbase --origfmt SRR3614678
fastq-dump --gzip --dumpbase --origfmt SRR3614679
fastq-dump --gzip --dumpbase --origfmt SRR3614680

## SE RNA barcodes, rep 1 and 2
fastq-dump --gzip --dumpbase --origfmt SRR3614681
fastq-dump --gzip --dumpbase --origfmt SRR3614682
fastq-dump --gzip --dumpbase --origfmt SRR3614683
fastq-dump --gzip --dumpbase --origfmt SRR3614684

