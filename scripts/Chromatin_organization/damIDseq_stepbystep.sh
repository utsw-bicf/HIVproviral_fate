#! /bin/bash
### Process DamID-seq
### Following https://github.com/sarahhcarl/flychip/wiki/Basic-DamID-analysis-pipeline

#SBATCH --job-name=damID_sbs
#SBATCH --partition=super
#SBATCH --output=sbs.%j.out
#SBATCH --error=sbs.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

fastqs="SRR5261759 SRR5261760 SRR5261761 SRR5261762"

###Trim DamID adapters
echo "#################### Trim ####################"
module load cutadapt/1.9.1
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/TrimReads

for fastq in ${fastqs}; do
  cutadapt -a GATCCTCGGCCGCGACC \
    -g ^GGTCGCGGCCGAGGATC \
    -o /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/TrimReads/${fastq}.fq.gz \
    /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/${fastq}.fastq.gz
done

### Assess read qualities
echo "#################### Quality Check ####################"
module load fastqc/0.11.5
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/fastqc

for fastq in ${fastqs}; do
  fastqc /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/TrimReads/${fastq}.fq.gz \
    --outdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/fastqc
done

### Map reads with bowtie2
echo "#################### bowtie2 ####################"
module load bowtie2/2.3.2
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped

for fastq in ${fastqs}; do
  bowtie2 -x GRCh38 \
    -U /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/TrimReads/${fastq}.fq.gz \
    -S /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastq}.sam
done

### Convert to bam, sort, and index reads
echo "#################### Samtools ####################"
module load samtools

for fastq in ${fastqs}; do
 samtools view -bS /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}.sam 
    -o /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}_unsorted.bam  
 samtools sort /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}_unsorted.bam ${fastqs}  
 samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}.bam  
 samtools idxstats /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}.bam
done

### Convert to bed, extend reads
echo "#################### bedtools ####################"
module load bedtools/2.26.0

for fastq in ${fastqs}; do
  bamToBed -i /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}.bam \
    >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}.bed
  sort -k1,1 /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}.bed \
    >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}.sorted.bed
  bedtools slop -l 150 \
    -r 0 \
    -s \
    -i /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}.sorted.bed \
    -g /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes \
    >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}.ext200.bed
done

### Visualize reads
echo "#################### bw ####################"
module load UCSC_userApps/v317

for fastq in ${fastqs}; do
  bedtools genomecov -i /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}.ext200.bed \
    -g /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes -d | awk -vOFS='\t' '{print $1,$2,$2+1,$3}' | gzip >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}.ext200.density.gz
  gunzip -c /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}.ext200.density.gz | awk -vOFS='\t' '($4!=0){if(!chrom[$1]){print "variableStep chrom="$1;chrom[$1]=1};print $2,$4}' | gzip >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}.ext200.wig.gz
  wigToBigWig /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}.ext200.wig.gz /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}.ext200.bw
done


### Find GATC sites
echo "#################### HOMER ####################"
module load homer/4.9

seq2profile.pl GATC 0 GATC >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/output.motif
scanMotifGenomeWide.pl GATC /project/shared/bicf_workflow_ref/human/GRCh38/genome.fa -bed | awk 'BEGIN {OFS="\t"}; {print $1, $8, $9, $10}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/GATC_sites.bed

### Make tables of read counts in all GATC fragments
echo "#################### counts ####################"
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/DE

for fastq in ${fastqs}; do
  bedtools coverage \
    -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/mapped/${fastqs}.ext200.bed \
    -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/GATC_sites.bed \
    > /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/DE/${fastqs}_GATC_counts.bed
done

echo "#################### END OF SCRIPT ####################"
