#! /bin/bash

#SBATCH --job-name=dastk_loop
#SBATCH --partition=super
#SBATCH --output=dastk_loop.%j.out
#SBATCH --error=dastk_loop.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load python/2.7.x-anaconda

### Get the left end of loop
awk '{OFS="\t"}; NR > 1 && $1==$4 {print $1, $2, $3}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/left_loop.bed

### Get the right end of loop
awk '{OFS="\t"}; NR > 1 && $1==$4 {print $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/right_loop.bed


### Run dastk (git clone git@git.biohpc.swmed.edu:BICF/Astrocyte/DAStk_modified.git)
### Jaspar motif files are located: /project/shared/bicf_workflow_ref/dastk_JASPAR_motifs/reformatted_files_human
python /work/BICF/s185797/programs/dastk/DAStk_modified/DAStk/process_atac.py \
  -x left_loop_jaspar \
  -e /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/left_loop.bed \
  -m /project/shared/bicf_workflow_ref/dastk_JASPAR_motifs/reformatted_files_human \
  -t 8

python /work/BICF/s185797/programs/dastk/DAStk_modified/DAStk/process_atac.py \
  -x right_loop_jaspar \
  -e /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/right_loop.bed \
  -m /project/shared/bicf_workflow_ref/dastk_JASPAR_motifs/reformatted_files_human \
  -t 8

### HocoMoco
python /work/BICF/s185797/programs/dastk/DAStk_modified/DAStk/process_atac2.py \
  -x left_loop_HocoMoco \
  -e /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/left_loop.bed \
  -m /project/shared/bicf_workflow_ref/dastk_JASPAR_motifs/hocomoco/HOCOMOCO_v11_p1e-6_grch38 \
  -t 8

python /work/BICF/s185797/programs/dastk/DAStk_modified/DAStk/process_atac2.py \
  -x right_loop_HocoMoco \
  -e /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/right_loop.bed \
  -m /project/shared/bicf_workflow_ref/dastk_JASPAR_motifs/hocomoco/HOCOMOCO_v11_p1e-6_grch38 \
  -t 8


### Sort from low to high scores
awk -F "," '{OFS="\t"} {print $1, $2}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/left_loop_jaspar_md_scores.txt | sort -rn -k 2,2 >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/left_loop_jaspar_md_scores_sorted.txt

awk -F "," '{OFS="\t"} {print $1, $2}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/right_loop_jaspar_md_scores.txt | sort -rn -k 2,2 >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/right_loop_jaspar_md_scores_sorted.txt

awk -F "," '{OFS="\t"} {print $1, $2}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/left_loop_HocoMoco_md_scores.txt | sort -rn -k 2,2 >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/left_loop_HocoMoco_md_scores_sorted.txt

awk -F "," '{OFS="\t"} {print $1, $2}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/right_loop_HocoMoco_md_scores.txt | sort -rn -k 2,2 >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/dastk/right_loop_HocoMoco_md_scores_sorted.txt
