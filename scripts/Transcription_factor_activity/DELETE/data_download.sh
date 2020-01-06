#! /bin/sh
### Download HiChiP data for d'Orso lab

#SBATCH --job-name=hiChiP
#SBATCH --partition=super
#SBATCH --output=hiChiP.%j.out
#SBATCH --error=hiChiP.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load sra_toolkit/2.8.2-1

### HiChiP, PE
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR601/SRR6010260/SRR6010260.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR601/SRR6010261/SRR6010261.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/raw_data/
