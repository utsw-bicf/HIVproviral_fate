#! /bin/sh

########## This script uses already done subcompartments from Hi-C data
### From "Rao SS, Huntley MH, Durand NC, Stamenova EK, Bochkov ID, Robinson JT, Sanborn AL, Machol I, Omer AD, Lander ES, et al. 2014. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. Cell 159: 1665–1680."

########## Download the subcompartments
#ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_subcompartments.bed.Egz
zcat GSE63525_GM12878_subcompartments.bed.gz >GSE63525_GM12878_subcompartments.bed

########## Genome used was b37; use lift over to hg38
### hg19 is the same as b37 except for the nameing of chromosomes (b37:1, hg19:chr1
### directly
### chain file: https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/funcotator/data_sources/gnomAD/b37ToHg38.over.chain
### chain file: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
liftOver GSE63525_GM12878_subcompartments.bed hg19ToHg38.over.chain.gz GSE63525_GM12878_subcompartments_hg38.bed GSE63525_GM12878_subcompartments_hg38_notlifted.bed

########## Do chromHMM with damID seq using this lift over
awk '{OFS="\t"};$4=="A1"{print $1,$2,$3}' GSE63525_GM12878_subcompartments_hg38.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM/GM12878_A1.bed
awk '{OFS="\t"};$4=="A2"{print $1,$2,$3}' GSE63525_GM12878_subcompartments_hg38.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM/GM12878_A2.bed
awk '{OFS="\t"};$4=="B1"{print $1,$2,$3}' GSE63525_GM12878_subcompartments_hg38.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM/GM12878_B1.bed
awk '{OFS="\t"};$4=="B2"{print $1,$2,$3}' GSE63525_GM12878_subcompartments_hg38.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM/GM12878_B2.bed
awk '{OFS="\t"};$4=="B3"{print $1,$2,$3}' GSE63525_GM12878_subcompartments_hg38.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM/GM12878_B3.bed

# Run chromHMM OverlapEnrichment
### Re-run with enhancer and superenhancer databases
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_ases_dam/OverlapEnrichment_enhancers_subcompartments
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_ases_dam/state15/jurkat_15_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone_ases_dam/OverlapEnrichment_enhancers_subcompartments/overlap_enhancers_subcompartments
