#! /bin/sh
### Download Bhive data for d'Orso lab

#SBATCH --job-name=dd
#SBATCH --partition=super
#SBATCH --output=dd.%j.out
#SBATCH --error=dd.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load sra_toolkit/2.8.2-1

### Bhive
### PE Integration map, rep 1 and 2
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR361/SRR3614675/SRR3614675.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR361/SRR3614676/SRR3614676.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/

## SE DNA barcodes, rep 1 and 2
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR361/SRR3614677/SRR3614677.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR361/SRR3614678/SRR3614678.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR361/SRR3614679/SRR3614679.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR361/SRR3614680/SRR3614680.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/

## SE RNA barcodes, rep 1 and 2
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR361/SRR3614681/SRR3614681.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR361/SRR3614682/SRR3614682.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR361/SRR3614683/SRR3614683.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR361/SRR3614684/SRR3614684.sra -P /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/

