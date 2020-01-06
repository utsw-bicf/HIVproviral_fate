#! /bin/bash

#SBATCH --job-name=cedbHOMER
#SBATCH --partition=super
#SBATCH --output=cedbHOMER.%j.out
#SBATCH --error=cedbHOMER.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load homer/4.9
module load samtools
module load R

### only use regular chromosomes
loc='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/filterReads'
for item in SRR2043614 SRR1603650 SRR1603654 GSM1603211 SRR2043612 SRR1603649 SRR1603652 GSM1603229; do
  samtools view -b ${loc}/${item}.filt.nodup.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM >${loc}/${item}.filt.nodup_chronly.bam
  samtools index ${loc}/${item}.filt.nodup_chronly.bam
done

### Create tag directory for H3K27ac
makeTagDirectory /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27_HOMER \
  ${loc}/SRR2043614.filt.nodup_chronly.bam \
  ${loc}/SRR1603650.filt.nodup_chronly.bam \
  ${loc}/SRR1603654.filt.nodup_chronly.bam \
  ${loc}/GSM1603211.filt.nodup_chronly.bam \
  -format sam  

### Create an input tag directory for input of H3K27ac
makeTagDirectory /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/INPUT_H3K27_HOMER \
  ${loc}/SRR2043612.filt.nodup_chronly.bam \
  ${loc}/SRR1603649.filt.nodup_chronly.bam \
  ${loc}/SRR1603652.filt.nodup_chronly.bam \
  ${loc}/GSM1603229.filt.nodup_chronly.bam \
  -format sam  
  
### Find enhancers

### Find superenhancers
findPeaks /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27_HOMER \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/INPUT_H3K27_HOMER \
  -style super \
  -typical H3K27_HOMER/enhancers.txt \
  -o auto

### Separate enhancer database into enhancers and promoters
### Easiest way is to annotate the peaks using chipseaker
grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27_HOMER/enhancers.txt | awk '{OFS="\t"}; {print $2,$3,$4,$1,$8,$5,$6,$7,$9,$11}' >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27_HOMER/enhancers.narrowpeak
echo -e "Condition\tPeaks
enhancers\t/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27_HOMER/enhancers.narrowpeak" >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27_HOMER/design_annotate_enhancers.tsv

Rscript /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/scripts/annotate_peaks.R /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27_HOMER/design_annotate_enhancers.tsv GRCh38

# Fix file
perl -pi -e "s|3\' UTR|3\'UTR|g;s|5\' UTR|5\'UTR|g;s|Distal Intergenic|DistalIntergenic|g;s|Downstream \(.*\)|Downstream|g;s|Promoter \(1-2kb\)|Promoter\(1-2kb\)|g;s|Promoter \(<=1kb\)|Promoter\(<=1kb\)|g;s|Promoter \(2-3kb\)|Promoter\(2-3kb\)|g;s|Intron \(.*\)|Intron|g;s|Exon \(.*\)|Exon|g" /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27_HOMER/enhancers.chipseeker_annotation.tsv

# Promoter <= 2kb from TSS (column 20)
awk -F "\t" '$20 < 2000 && $20 > -2000 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27_HOMER/enhancers.chipseeker_annotation.tsv >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27_HOMER/HOMER_promoter_2kb.tsv

# enhancer
awk  -F "\t" '$20 > 2000 || $20 < -2000 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27_HOMER/enhancers.chipseeker_annotation.tsv | sed '1d' >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27_HOMER/HOMER_enhancer.tsv

### Run chromHMM with refseq, enhancers and superenhancers
# Transfer over refseq stuff from chromHMM
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM
ln -s /work/BICF/s185797/programs/ChromHMM/COORDS/hg38/*.bed.gz /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM

# make promoter, enhancer, and superenhancer bed file (chr, start,end)
awk '{OFS="\t"}; {print $1, $2, $3}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27_HOMER/HOMER_promoter_2kb.tsv >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM/H_promoter.bed
awk '{OFS="\t"}; {print $1, $2, $3}' /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27_HOMER/HOMER_enhancer.tsv >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM/H_enhancers.bed
grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27_HOMER/superEnhancers.txt | awk '{OFS="\t"}; {print $2, $3, $4}' >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM/H_superenhancers.bed

# Run chromHMM OverlapEnrichment
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/OverlapEnrichment_enhancers
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/OverlapEnrichment_enhancers/overlap_enhancers
