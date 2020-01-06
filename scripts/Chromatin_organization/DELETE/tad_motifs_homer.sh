#!/bin/bash

#SBATCH --job-name=motif
#SBATCH --partition=super
#SBATCH --output=motif.%j.out
#SBATCH --error=motif.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

########## This script uses homer to identify motifs at the end of loops
module load homer/4.9

### Fix tads
## Left
grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.tad.2D.bed | awk '{OFS = "\t"}{print $1, $2, $2+15000, $7, $8, $9}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/left_15000_tad.bed

## Right
grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.tad.2D.bed | awk '{OFS = "\t"}{print $1, $3-15000, $3, $7, $8, $9}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/right_15000_tad.bed

### Set the genome
cp /project/shared/bicf_workflow_ref/human/GRCh38/genome.fa /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/motifs_homer/hg38.fa

findMotifsGenome.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/left_15000_tad.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/motifs_homer/hg38.fa \
  tad_Left_15000_homer_motif_output \
  -size given \
  -mask

findMotifsGenome.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/right_15000_tad.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/motifs_homer/hg38.fa \
  tad_Right_15000_homer_motif_output \
  -size given \
  -mask


######################################################################################################
## Left
grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.tad.2D.bed | awk '{OFS = "\t"}{print $1, $2, $2+8000, $7, $8, $9}' | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/left_8000_tad.bed
findMotifsGenome.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/left_8000_tad.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/motifs_homer/hg38.fa \
  tad_Left_8000_homer_motif_output \
  -size given \
  -mask

