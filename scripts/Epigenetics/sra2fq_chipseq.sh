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

#SRAs=`ls *.sra`
SRAs="SRR3722566"

for sra in $SRAs; do
fastq-dump --gzip --dumpbase --origfmt ${sra}
done
