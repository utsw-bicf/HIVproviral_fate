#! /bin/sh
### Run Hi-C data through HOMER

#SBATCH --job-name=HOMER
#SBATCH --mem 126976
#SBATCH --output=HOMER.%j.out
#SBATCH --error=HOMER.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### Follow: http://homer.ucsd.edu/homer/interactions2/quickndirty.html
module load homer/4.9

### Concatenate fastqs by library
#echo "Start concatenating libraries"
## lib1
#cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8244644_1.fastq.gz \
#/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8244645_1.fastq.gz \
#>/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489136_1.fastq.gz

#cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8244644_2.fastq.gz \
#/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8244645_2.fastq.gz \
#>/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489136_2.fastq.gz

## lib2
#cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8244646_1.fastq.gz \
#/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388751_1.fastq.gz \
#/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388752_1.fastq.gz \
#/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388753_1.fastq.gz \
#/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388754_1.fastq.gz \
#/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388755_1.fastq.gz \
#/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388756_1.fastq.gz \
#/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388757_1.fastq.gz \
#>/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489137_1.fastq.gz

#cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8244646_2.fastq.gz \
#/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388751_2.fastq.gz \
#/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388752_2.fastq.gz \
#/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388753_2.fastq.gz \
#/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388754_2.fastq.gz \
#/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388755_2.fastq.gz \
#/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388756_2.fastq.gz \
#/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388757_2.fastq.gz \
#>/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489137_2.fastq.gz
#echo "End concatenating libraries"

### Trim reads
#echo "Start read trimming"
#fastqs='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489136_1.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489136_2.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489137_1.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489137_2.fastq.gz'

#for fastq in ${fastqs}; do
#  homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 ${fastq}
#done
#echo "End read trimming"

### Make reference genome for bowtie using only canonical chromosomes (1-22,X,Y,M)
#echo "Start making canonical reference"
#mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/canonical_genome

#module load samtools
#samtools faidx /project/shared/bicf_workflow_ref/human/GRCh38/genome.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/canonical_genome/canonical_genome.fa
#module rm samtools

#module load bowtie2/2.2.8-intel
#mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/canonical_genome/bowtie2_index/
#bowtie2-build -f /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/canonical_genome/canonical_genome.fa /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/canonical_genome/bowtie2_index/genome
#echo "End making canonical reference"

### Alignment with bowtie2
#echo "Start alignment"
#mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/alignments
#tfastqs='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489136_1.fastq.gz.trimmed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489136_2.fastq.gz.trimmed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489137_1.fastq.gz.trimmed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489137_2.fastq.gz.trimmed'

#for tfastq in ${tfastqs}; do
#  outname=$(basename ${tfastq} .fastq.gz.trimmed)
#  bowtie2 -p 20 -x /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/canonical_genome/bowtie2_index/genome -U ${tfastq} >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/alignments/${outname}.sam
#done
#module rm bowtie2/2.2.8-intel
#echo "End alignment"

### HI-C tag directory
echo "Start Hi-C tagging"
#mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/HicHomerTagDir/
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/HicHomerTagDir/lib1
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/HicHomerTagDir/lib2
makeTagDirectory HicHomerTagDir/lib1/ alignments/GSM3489136_1.sam,alignments/GSM3489136_2.sam -tbp 1 -genome hg38 -checkGC -restrictionSite GATC -removePEbg -removeSpikes 10000 5
makeTagDirectory HicHomerTagDir/lib2 alignments/GSM3489137_1.sam,alignments/GSM3489137_2.sam -tbp 1 -genome hg38 -checkGC -restrictionSite GATC -removePEbg -removeSpikes 10000 5
echo "End Hi-C tagging"



