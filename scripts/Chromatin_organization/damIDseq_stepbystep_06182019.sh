#! /bin/bash
### Process DamID-seq
### Following "Dam mutants provide improved sensitivity and spatial resolution for profiling transcription factor binding"

#SBATCH --job-name=damID_sbs
#SBATCH --partition=super
#SBATCH --output=sbs.%j.out
#SBATCH --error=sbs.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### Find GATC sites
#echo "#################### HOMER ####################"
#module load homer/4.9

#seq2profile.pl GATC 0 GATC >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/output.motif
#perl -p -e 's|>(\w+).*$|>$1|g' /project/shared/bicf_workflow_ref/human/GRCh38/genome.fa >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/test.fa
#scanMotifGenomeWide.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/output.motif /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/test.fa -bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/GATC_sites.bed
#rm /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/test.fa

fastqs="SRR5261759 SRR5261760 SRR5261761 SRR5261762"

###Trim DamID adapters
#echo "#################### Trim ####################"
#module load cutadapt/1.9.1
#mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/TrimReads

#for fastq in ${fastqs}; do
#  cutadapt -a GATCCTCGGCCGCGACC \
#    -g ^GGTCGCGGCCGAGGATC \
#    -o /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/TrimReads/${fastq}.fq.gz \
#    /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/${fastq}.fastq.gz
#done

### Assess read qualities
#echo "#################### Quality Check ####################"
#module load fastqc/0.11.5
#mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/fastqc

#for fastq in ${fastqs}; do
#  fastqc /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/TrimReads/${fastq}.fq.gz \
#    --outdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/fastqc
#done

### Map reads with bowtie2
#echo "#################### bowtie2 ####################"
#module load bowtie2/2.3.2
#mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped

#for fastq in ${fastqs}; do
#  bowtie2 -x GRCh38 \
#    -U /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/TrimReads/${fastq}.fq.gz \
#    -S /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}.sam
#done

### Convert to bam, sort, keep only chr1-22,X,Y and index reads, remove unwanted files
#echo "#################### Samtools ####################"
#module load samtools

#for fastq in ${fastqs}; do
# samtools view -bS /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}.sam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}_unsorted.bam  
# samtools sort /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}_unsorted.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}_sorted.bam
# samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}_sorted.bam  
# samtools view -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}_sorted.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
# samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}.bam
# samtools idxstats /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}.bam >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}SRR5261759.idxstats.txt
#done
#rm /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/*.sam* /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/*_unsorted.bam* /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/*_sorted.bam*

### Convert to bed keeping, compare to GATC sites(perl script), reduce bam files
#echo "#################### bedtools ####################"
#module load bedtools/2.26.0
#module load picard/2.10.3
#module load samtools/intel/1.3

#mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/counts

## Make GATC file to intersect for multimapped reads
#awk '{OFS="\t"};{print $1,$2-2,$3+1,$4,$5,$6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/GATC_sites.bed | awk '{OFS="\t"}; $2 >1 {print $0}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/GATC_intersect.bed

#for fastq in ${fastqs}; do
  #bamToBed -i /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}.bam | sort -k1,1 -k 2,2n /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}.bed
  #perl /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Chromatin_organization/damIDseq_keepreads.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/GATC_sites.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}.keeplist.txt
#  java -jar /cm/shared/apps/picard/2.10.3/picard.jar FilterSamReads I=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}.bam O=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}_delete.bam READ_LIST_FILE=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}.keeplist.txt FILTER=includeReadList
#  samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}_delete.bam
#  bedtools intersect -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/${fastq}_delete.bam \
#    -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/GATC_intersect.bed \
#    -wa >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/counts/${fastq}_filt.bam
 # samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/counts/${fastq}_filt.bam
#done

### call reads as peaks
module load macs/2.1.0-20151222
#echo "#################### START calling peaks ####################"

macs2 callpeak \
  -t /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/counts/SRR5261760_filt.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/counts/SRR5261762_filt.bam \
  -c /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/counts/SRR5261759_filt.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/counts/SRR5261761_filt.bam \
  -n damIDseq \
  -g hs \
  --keep-dup all \
  --bw 300 \
  --qvalue 0.05 \
  --mfold 5 50 \
  --broad \
  --broad-cutoff 0.1 \
  -B

### Filter for significant
awk '{OFS="\t"}; $9>1.3 && $7 > 2 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/call_peaks/damIDseq_peaks.broadPeak >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/call_peaks/damIDseq_peaks_sig.broadPeak

### Convert bedgraph to bw
module load UCSC_userApps/v317 

macs2 bdgcmp -t damIDseq_treat_pileup.bdg \
  -c damIDseq_control_lambda.bdg \
  -o damIDseq_FE.bdg \
  -m FE

bedSort /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/call_peaks/damIDseq_FE.bdg /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/call_peaks/damIDseq_sorted.bedGraph

bedGraphToBigWig /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/call_peaks/damIDseq_sorted.bedGraph /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/call_peaks/damIDseq_signal.bw

#echo "#################### END calling peaks ####################"


########## Don't use below
### call peaks on unfiltered reads
#echo "#################### START calling peaks ####################"

#macs2 callpeak \
#  -t /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/SRR5261760.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/SRR5261762.bam \
#  -c /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/SRR5261759.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/mapped/SRR5261761.bam \
#  -n damIDseq_unfilt \
#  -g hs \
#  --keep-dup all \
#  --bw 300 \
#  --qvalue 0.05 \
#  --mfold 5 50 \
#  --broad \
#  --broad-cutoff 0.1

### Filter for significant
#awk '{OFS="\t"}; $9>1.3 && $7 > 2 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/call_peaks/damIDseq_unfilt_peaks.broadPeak >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/call_peaks/damIDseq_unfilt_peaks_sig.broadPeak
#echo "#################### END calling peaks ####################"
