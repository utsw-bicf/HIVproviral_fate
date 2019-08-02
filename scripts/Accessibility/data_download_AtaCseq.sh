#! /bin/sh
### Download RNAseq data for d'Orso lab
### Chipseq

#SBATCH --job-name=dl_Atac
#SBATCH --partition=super
#SBATCH --output=dla.%j.out
#SBATCH --error=dla.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load sra_toolkit/2.8.2-1

SRAs="SRR5720369 SRR5720371"

for sra in $SRAs; do
  sixsra=`echo ${sra} | cut -c1-6`
  wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${sixsra}/${sra}/${sra}.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/raw_data

  echo ${sra}.sra >>atac_checksum.txt
  vdb-validate ${sra}.sra >>atac_checksum.txt

  fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip --origfmt ${sra}
  rm /home2/s185797/ncbi/public/sra/${sra}.sra
done

