#! /bin/bash

#SBATCH --job-name=cedb
#SBATCH --partition=super
#SBATCH --output=cedb.%j.out
#SBATCH --error=cedb.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load samtools
module load subread/1.6.1 stringtie/1.3.2d-intel

##### This script will create an enhancer database from the H3K27ac mark
### Steps:
### 1. Merge bam files
### 2. separate chip peaks into promoter and enhancers
### 3. Run Rose (http://younglab.wi.mit.edu/super_enhancer_code.html) on enhancers only 
#####

### Step 1: Merge H3K27ac bam files
echo "Start Step 1: Merge H3K27ac bam files"
loc='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/filterReads'
 
samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_merged.bam ${loc}/SRR2043614.filt.nodup.bam ${loc}/SRR1603650.filt.nodup.bam ${loc}/SRR1603654.filt.nodup.bam ${loc}/GSM1603211.filt.nodup.bam

samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_merged.bam

samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Input_H3K27ac_merged.bam ${loc}/SRR2043612.filt.nodup.bam ${loc}/SRR1603649.filt.nodup.bam ${loc}/SRR1603652.filt.nodup.bam ${loc}/GSM1603229.filt.nodup.bam

samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Input_H3K27ac_merged.bam

echo "End Step 1: Merge H3K27ac bam files"

### Step 2: Separate H3K27ac consensus, replicated narrow peaks and separate into promoter (2kb +/- TSS) or other
### Take "other" as a "feature" and run featureCounts on the TTseq 
### Plot histogram of fpkm of TTseq (expect bimodal); use a cutoff to identify enhancers
### Note: I still keep lincRNA as an enhancer

echo "Start Step2: make enhancer database"

# Fix file
perl -p -e "s|3\' UTR|3\'UTR|g;s|5\' UTR|5\'UTR|g;s|Distal Intergenic|DistalIntergenic|g;s|Downstream \(.*\)|Downstream|g;s|Promoter \(1-2kb\)|Promoter\(1-2kb\)|g;s|Promoter \(<=1kb\)|Promoter\(<=1kb\)|g;s|Promoter \(2-3kb\)|Promoter\(2-3kb\)|g;s|Intron \(.*\)|Intron|g;s|Exon \(.*\)|Exon|g" /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/peakAnnotation/H3K27ac.chipseeker_annotation.tsv >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac.chipseeker_annotation_edited.tsv

# Promoter <= 2kb from TSS (column 20)
awk -F "\t" '$20 < 2000 && $20 > -2000 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac.chipseeker_annotation_edited.tsv >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_promoter_2kb.chipseeker_annotation.tsv

# Other (not promoter)
awk  -F "\t" '$20 > 2000 || $20 < -2000 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac.chipseeker_annotation_edited.tsv | perl -p -e "s|Promoter\(2-3kb\)|farPromoter|g;s|3\'|3prime|g;s|5\'|5prime|g" >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_not_promoter.chipseeker_annotation.tsv

# Make gtf file for featureCounts
awk  -F "\t" '{OFS="\t"}; FNR>1 {print $1, "chipseeker", "notpromoter", $2, $3, ".", ".", ".", "gene_id \""$5"\"; gene_type \""$12"\"; gene_status \"KNOWN\"; gene_name \""$21"\";"}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_not_promoter.chipseeker_annotation.tsv >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_not_promoter.gtf

sed -i '1i\##description: not evidence-based notpromoters of the human genome (GRCh38)\n##provider: GENCODE\n##contact: holly.ruess@utsouthwestern.edu\n##format: gtf\n##date: 2019-05-17' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_not_promoter.gtf

# Combine TTseq bam files
samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Ttseq_CRAMER_merged.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.dedup.bam

samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Ttseq_CRAMER_merged.bam

# Run featureCounts on TTseq
featureCounts -s 0 -T 4 -p -g gene_id -t notpromoter -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_not_promoter.gtf -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Ttseq_CRAMER.cts /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Ttseq_CRAMER_merged.bam

# Run strintie
# Don't use:stringtie /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Ttseq_CRAMER_merged.bam -G /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_not_promoter.gtf -B -e -o denovo.gtf -A /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Ttseq_CRAMER.fpkm.txt

# Take the string tie output and identify peaks that have a count > 200
awk '$7 > 200 {print $1}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Ttseq_CRAMER.cts >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Ttseq_peaks.txt

# Take significant, replicated enchancers and overlap with ttseq_peaks
grep -wFf /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Ttseq_peaks.txt /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_not_promoter.chipseeker_annotation.tsv | awk '{OFS="\t"}; $10 > 1.30102999566 {print $0}' >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_ttSeq_possible_enhancers.tsv

# Make file gtf for rose
# Rose requirements:
	# .gff must have the following columns:
	# 1: chromosome (chr#)
	# 2: unique ID for each constituent enhancer region
	# 4: start of constituent
	# 5: end of constituent
	# 7: strand (+,-,.)
	# 9: unique ID for each constituent enhancer region
	# NOTE: if value for column 2 and 9 differ, value in column 2 will be used
awk '{OFS="\t"}; {print $1,$5,$12,$2,$3,$6,$7,$7,$5}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_ttSeq_possible_enhancers.tsv >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_ttSeq_possible_enhancers.gff

echo "End Step2: make enhancer database"

### 3. Run Rose (http://younglab.wi.mit.edu/super_enhancer_code.html) on enhancers only 
### Install Rose: module load python/2.7.x-anaconda; pip install --user pyrose
### nano /home2/s185797/.local/bin/pyrose; and comment out import ipsea
#python ROSE_main.py -g GENOME_BUILD -i INPUT_CONSTITUENT_GFF -r RANKING_BAM -o OUTPUT_DIRECTORY
pyrose -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_ttSeq_possible_enhancers.gff \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_merged.bam \
  -c /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Input_H3K27ac_merged.bam \
  -a HG38 \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/rose_enhancer_results


### Install Rose: git clone https://bitbucket.org/young_computation/rose.git
cd /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/rose/rose
python ROSE_main.py \
  -g hg38 \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_ttSeq_possible_enhancers.gff \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_merged.bam \
  -c /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Input_H3K27ac_merged.bam \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/rose_enhancer_results \
  -s 12500 \
  -t 2000
