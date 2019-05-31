#! /bin/bash

#SBATCH --job-name=cedb_rose
#SBATCH --partition=super
#SBATCH --output=cedb_rose.%j.out
#SBATCH --error=cedb_rose.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load samtools
module load python/2.7.6-epd

python ROSE_main.py \
  -g HG18 \
  -i /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_ttSeq_possible_enhancers.gff \
  -r /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/H3K27ac_merged.bam \
  -c /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Input_H3K27ac_merged.bam \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/rose_enhancer_results \
  -s 12500 \
  -t 2500

### Filter for superenhancer vs enhancer
# Superenhancer
grep -v "#\|REGION_ID" /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/rose_enhancer_results/H3K27ac_ttSeq_possible_enhancers_AllEnhancers.table.txt | awk '{OFS="\t"}; $10 ~ /0/ {print $2, $3, $4}' >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM/R_superenhancers.bed
# Enhancers
grep -v "#\|REGION_ID" /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/rose_enhancer_results/H3K27ac_ttSeq_possible_enhancers_AllEnhancers.table.txt | awk '{OFS="\t"}; $10 ~ /1/ {print $2, $3, $4}' >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/database_for_chromHMM/R_enhancers.bed
