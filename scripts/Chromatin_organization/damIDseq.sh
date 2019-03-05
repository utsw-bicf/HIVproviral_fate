#! /bin/sh
### Process DamID-seq
### Used already made pipeline https://owenjm.github.io/damidseq_pipeline/ v 1.4.5 downloaded 2/25/19
### Reference: Marshall OJ and Brand AH. (2015) damidseq_pipeline: an automated pipeline for processing DamID sequencing datasets. Bioinformatics. Oct 15;31(20):3371-3. doi: 10.1093/bioinformatics/btv386

#SBATCH --job-name=damID
#SBATCH --partition=super
#SBATCH --output=damID.%j.out
#SBATCH --error=damID.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load bowtie2/2.3.2
module load samtools/1.6
#module load perl/5.28.0
module load bedtools/2.26.0

### Change the first line of the damidseq_pipeline to look for perl in a general way

### Build your own GATC file
#perl /work/BICF/s185797/programs/owenjm-damidseq_pipeline-4ec63e9/gatc.track.maker.pl --name=GRCh38_GATC /project/shared/bicf_workflow_ref/GRCh38/genome.fa

### Permanatey set the pipeline options; not sure if this worked, so specify it on the real command
#perl /work/BICF/s185797/programs/owenjm-damidseq_pipeline-4ec63e9/damidseq_pipeline --gatc_frag_file=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/GRCh38_GATC.GATC.gff \
#  --bowtie2_genome_dir=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/bowtie2_index/GRCh38

### Make softlink for the program
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261759.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261759.Index12.fastq.gz
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261760.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261760.Index5.fastq.gz

### Make index file for program
#echo -e "A5\tpolII" >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/index.txt
#echo -e "A12\tDam" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/index.txt

### Run program
#perl /work/BICF/s185797/programs/owenjm-damidseq_pipeline-4ec63e9/damidseq_pipeline --gatc_frag_file=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/GRCh38_GATC.GATC.gff \
#  --bowtie2_genome_dir=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/bowtie2_index/GRCh38 \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261759.Index12.fastq.gz \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261760.Index5.fastq.gz

### Results from above are in folder and will be deleted
############################################################################################################
############################################################################################################
#####START AGAIN
#SRR5261759, ctl
#SRR5261760, LaminB1
#SRR5261761, ctl
#SRR5261762, LaminB1

### Changes to original script
##Line 1: #!/usr/bin/env perl -w

### Reset any previous defaults
perl /work/BICF/s185797/programs/owenjm-damidseq_pipeline-4ec63e9/damidseq_pipeline --reset_defaults 

### Make softlink for the program
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261759.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261759.Index1.fastq.gz
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261760.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261760.Index2.fastq.gz
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261761.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261761.Index3.fastq.gz
#ln -s /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261762.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261762.Index4.fastq.gz

### Make index file for program
echo -e "A1\tDam1" >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/index.txt
echo -e "A2\tLaminB11" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/index.txt


### Run program
perl /work/BICF/s185797/programs/owenjm-damidseq_pipeline-4ec63e9/damidseq_pipeline \
  --bins=500 \
  --norm_method=rpm \
  --gatc_frag_file=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/GRCh38_GATC.GATC.gff \
  --bowtie2_genome_dir=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/bowtie2_index/GRCh38 \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261759.Index1.fastq.gz \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261760.Index2.fastq.gz 


###Then run
#echo -e "A3\tDam2" >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/index.txt
#echo -e "A4\tLaminB12" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/index.txt

#perl /work/BICF/s185797/programs/owenjm-damidseq_pipeline-4ec63e9/damidseq_pipeline \
#  --bins=500 \
#  --norm_method=rpm \
#  --gatc_frag_file=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/GRCh38_GATC.GATC.gff \
#  --bowtie2_genome_dir=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/bowtie2_index/GRCh38 \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261761.Index3.fastq.gz \
 # /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR5261762.Index4.fastq.gz 
