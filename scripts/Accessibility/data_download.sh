#! /bin/sh
### Download DNAse and MNase seq data for d'Orso lab

#SBATCH --job-name=DM
#SBATCH --partition=super
#SBATCH --output=DM.%j.out
#SBATCH --error=DM.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### DNAse-seq
wget https://www.encodeproject.org/files/ENCFF001DPG/@@download/ENCFF001DPG.fastq.gz
wget https://www.encodeproject.org/files/ENCFF001DPF/@@download/ENCFF001DPF.fastq.gz

### MNase-seq, from Ivan
cp /archive/shared/DOrso_BICF/Jkt106_MNase_Seq_hg38/Raw_FASTQGZ/* .
