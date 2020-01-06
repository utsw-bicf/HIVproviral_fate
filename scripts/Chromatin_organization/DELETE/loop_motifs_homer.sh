#!/bin/bash

#SBATCH --job-name=motif
#SBATCH --partition=super
#SBATCH --output=motif.%j.out
#SBATCH --error=motif.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

########## This script uses homer to identify motifs at the end of loops
module load homer/4.9

### Set the genome
cp /project/shared/bicf_workflow_ref/human/GRCh38/genome.fa /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/motifs_homer/hg38.fa

findMotifsGenome.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/left_loop.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/motifs_homer/hg38.fa \
  Left_homer_motif_output \
  -size given \
  -mask

findMotifsGenome.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/right_loop.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/motifs_homer/hg38.fa \
  Right_homer_motif_output \
  -size given \
  -mask
