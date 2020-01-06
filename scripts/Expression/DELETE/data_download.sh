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
#wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR400/SRR4000396/SRR4000396.sra

### GSM2453329: SRR5170980
### single-end
#wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR517/SRR5170980/SRR5170980.sra

### 4sUseq and TTseq
###TTseq rep1-2
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR400/SRR4000388/SRR4000388.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR400/SRR4000389/SRR4000389.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/raw_data/

### 4sUseq, rep 1-3
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR213/SRR2130286/SRR2130286.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR213/SRR2130287/SRR2130287.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR213/SRR2130288/SRR2130288.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/raw_data/
