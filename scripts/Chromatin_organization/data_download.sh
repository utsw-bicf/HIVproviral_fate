#! /bin/sh
### Download Hi-C, DamID, and CTCF data for d'Orso lab

#SBATCH --job-name=CO
#SBATCH --partition=super
#SBATCH --output=CO.%j.out
#SBATCH --error=CO.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load sra_toolkit/2.8.2-1

### Hi-C, PE
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR824/SRR8244644/SRR8244644.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR824/SRR8244645/SRR8244645.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR824/SRR8244646/SRR8244646.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR838/SRR8388751/SRR8388751.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR838/SRR8388752/SRR8388752.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR838/SRR8388753/SRR8388753.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR838/SRR8388754/SRR8388754.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR838/SRR8388755/SRR8388755.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR838/SRR8388756/SRR8388756.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR838/SRR8388757/SRR8388757.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/

### DamID, SE
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR526/SRR5261759/SRR5261759.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR526/SRR5261760/SRR5261760.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR526/SRR5261761/SRR5261761.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR526/SRR5261762/SRR5261762.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/

### CTCF, SE
# Input
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR202/SRR2029842/SRR2029842.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/
# Experiment
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR202/SRR2029843/SRR2029843.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/

# Experiment, no input given
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR014/SRR014986/SRR014986.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/



























