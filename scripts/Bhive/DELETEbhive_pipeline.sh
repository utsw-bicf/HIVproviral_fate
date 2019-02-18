#! /bin/sh
### Step 1: barcodes

#SBATCH --job-name=rtg
#SBATCH --partition=super
#SBATCH --output=rtg.%j.out
#SBATCH --error=rtg.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### Load Modules
module load python/2.7.5
source ~/.bash_profile

### soft-link files into folder
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/SRR3614675_1.fastq.gz .
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/SRR3614675_2.fastq.gz .
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/SRR3614676_1.fastq.gz .
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/SRR3614676_2.fastq.gz .

### Extract barcodes
#python /work/BICF/s185797/programs/BHIVE_for_single_provirus_transcriptomics-master/src/mapping_extract_barcodes.py -b TATAGTGAGTCGTA -l AGCCCTTCCA -r CGCTTTTA --barcodes ipcr_1.bcd --reads ipcr_1.fastq SRR3614675_1.fastq.gz SRR3614675_2.fastq.gz &
#python /work/BICF/s185797/programs/BHIVE_for_single_provirus_transcriptomics-master/src/mapping_extract_barcodes.py -b TATAGTGAGTCGTA -l AGCCCTTCCA -r CGCTTTTA --barcodes ipcr_2.bcd --reads ipcr_2.fastq SRR3614676_1.fastq.gz SRR3614676_2.fastq.gz &
#wait

### Cluster barcodes
starcode -d 1 --seq-id -t 8 ipcr_1.bcd > ipcr_1.stc &
starcode -d 1 --seq-id -t 8 ipcr_2.bcd > ipcr_2.stc &
wait

python /work/BICF/s185797/programs/BHIVE_for_single_provirus_transcriptomics-master/src/mapping_barcode_seqnames.py ipcr_1.stc ipcr_1.fastq > ipcr_integ_site_1.fastq &
python /work/BICF/s185797/programs/BHIVE_for_single_provirus_transcriptomics-master/src/mapping_barcode_seqnames.py ipcr_2.stc ipcr_2.fastq > ipcr_integ_site_2.fastq &
wait
