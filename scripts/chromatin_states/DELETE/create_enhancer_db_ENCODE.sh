#! /bin/bash
### Use encode data to identify superenhancers

### CD4
#CL:0000897\|CL:0000624\|CL:0000792\|CL:0000905\|CL:0000895

### T-cell
#CL:0000084
#CL:0000899
module load bedtools/2.25.0

########## Parse the database from Chris
grep "CL:0000897\|CL:0000624\|CL:0000792\|CL:0000905\|CL:0000895" /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/get_encode_run_only_ChIP/GRCh38/metadata_GRCh38_20190523.tsv | grep "GRCh38" | grep -i "chip" | grep -i "peaks" | awk '{OFS="\t"} $2 != "bam" {print $0}' | grep -v "extremely low read depth" | awk -F "\t" '{OFS="\t"} $3 = "replicated peaks" {print $0}' | awk -F "\t" '{OFS="\t"} $31 ~ /2/ {print $0}' | grep -i "H3K27ac" >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/get_encode_run_only_ChIP/GRCh38/metadata_GRCh39_use.tsv

########## Make sliding window of reference genome
### Make chrom sizes file
grep -v "\_" /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/chrom.sizes

bedtools makewindows \
  -g /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/chrom.sizes \
  -w 12500 \
  -s 500 \
  >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/GRCh38_by12500_every500.bed

########## Make counts for each file from Chris's data
names="ENCFF625UHE ENCFF793UGX ENCFF868HON ENCFF747POK ENCFF689YGJ ENCFF083VQO ENCFF831PIG ENCFF952NOO ENCFF693ZID ENCFF348THZ ENCFF391ZDT ENCFF031NSV ENCFF522EGN ENCFF515AFW ENCFF802EZB ENCFF901AMN ENCFF831NID ENCFF857UAG ENCFF469LFE ENCFF355PHN ENCFF714CRK ENCFF338LAI"

for name in ${names}; do
bedtools intersect \
  -a /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/GRCh38_by12500_every500.bed \
  -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/get_encode_run_only_ChIP/GRCh38/ChIP-seq/H3K27ac/replicated_peaks/${name}.bed.gz \
  -c \
  >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/bed_counts/${name}.counts.bed
done

########## Put together the data in R
module load R
Rscript /home2/s185797/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/chromatin_states/create_enhancer_db_ENCODE.R
bedtools merge -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/bed_counts/unmerged_coordinates.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/bed_counts/GRCh38_ENCODE_superenhancers.bed

bedtools merge -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/bed_counts/unmerged_coordinates.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM/ENCODE_superenhancers.bed

########## run chromHMM
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/OverlapEnrichment_enhancers
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/learn_states/histone_learn15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone/OverlapEnrichment_enhancers/overlap_enhancers
