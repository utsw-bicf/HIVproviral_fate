#! /bin/bash

#SBATCH --job-name=dastk
#SBATCH --partition=super
#SBATCH --output=dastk.%j.out
#SBATCH --error=dastk.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load python/2.7.x-anaconda

### This file finds the MD score for motifs based on Hi and Low HIV expression
# Separated HIV expression for hi and low (0.5); extend bp by 50
awk '{OFS="\t"} $5 > 0.5 {print $1, $2-50, $3+50, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/HIV_SortedByExp.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/dastk_HIVexp/HIVexp_high.bed
awk '{OFS="\t"} $5 < -0.5 {print $1, $2-50, $3+50, $4, $5, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/HIV_SortedByExp.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/dastk_HIVexp/HIVexp_low.bed

### Run dastk (git clone git@git.biohpc.swmed.edu:BICF/Astrocyte/DAStk_modified.git)
### Jaspar motif files are located: /project/shared/bicf_workflow_ref/dastk_JASPAR_motifs/reformatted_files_human
python /work/BICF/s185797/programs/dastk/DAStk_modified/DAStk/process_atac.py \
  -x HIV_low \
  -e /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/dastk_HIVexp/HIVexp_low.bed \
  -m /project/shared/bicf_workflow_ref/dastk_JASPAR_motifs/reformatted_files_human \
  -t 8

python /work/BICF/s185797/programs/dastk/DAStk_modified/DAStk/process_atac.py \
  -x HIV_high \
  -e /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/dastk_HIVexp/HIVexp_high.bed \
  -m /project/shared/bicf_workflow_ref/dastk_JASPAR_motifs/reformatted_files_human \
  -t 8

### Need to make conda environment
#conda create -n dastk_env python=2.7

source activate dastk_env
#pip install adjustText
#pip install scipy
#pip install pandas


python /work/BICF/s185797/programs/dastk/DAStk_modified/DAStk/differential_md_score.py -x HIV -1 low -2 high -p 0.001 -b



########## HocoMoco
python /work/BICF/s185797/programs/dastk/DAStk_modified/DAStk/process_atac2.py \
  -x HM_low \
  -e /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/dastk_HIVexp/HIVexp_low.bed \
  -m /project/shared/bicf_workflow_ref/dastk_JASPAR_motifs/hocomoco/HOCOMOCO_v11_p1e-6_grch38 \
  -t 8

python /work/BICF/s185797/programs/dastk/DAStk_modified/DAStk/process_atac2.py \
  -x HM_high \
  -e /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/dastk_HIVexp/HIVexp_high.bed \
  -m /project/shared/bicf_workflow_ref/dastk_JASPAR_motifs/hocomoco/HOCOMOCO_v11_p1e-6_grch38 \
  -t 8

python /work/BICF/s185797/programs/dastk/DAStk_modified/DAStk/differential_md_score.py -x HM -1 low -2 high -p 0.005 -b





source deactivate
