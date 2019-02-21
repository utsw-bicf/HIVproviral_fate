#! /bin/sh
### Download HiChiP data for d'Orso lab

#SBATCH --job-name=HiChiP
#SBATCH --partition=super
#SBATCH --output=HiChiP.%j.out
#SBATCH --error=HiChiP.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load sra_toolkit/2.8.2-1

### HiChiP
### paired-end
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR6010260
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt SRR6010261


