#! /bin/bash

########## This script tries to identify the most common motifs in Loops
########## Tried CTCF but only explains a portion

module load bedtools

##### Combine left and right side loop into bed file
grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D.bed | awk '{OFS = "\t"} {print $1, $2, $3}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_motifs_common/left_loop.txt
grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D.bed | awk '{OFS = "\t"} {print $4, $5, $6}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_motifs_common/right_loop.txt
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_motifs_common/left_loop.txt /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_motifs_common/right_loop.txt | sort -k 1,1 -k 2,2n > /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_motifs_common/input_loop_leftright.txt


##### Loop through jaspar possibilities and calculate the line numbers
for jaspar in /project/shared/bicf_workflow_ref/dastk_JASPAR_motifs/reformatted_files_human/*.bed ; do
  name=$(basename -s .bed ${jaspar})
  count=$(bedtools intersect -wa -a /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_motifs_common/input_loop_leftright.txt -b <(grep -v "HLA" ${jaspar}) | wc -l)
  echo -e ${name}"\t"${count} >>/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_motifs_common/jaspar_motif_counts.txt
done

##### Sort file from high to low
sort -k 2,2nr /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/loop_motifs_common/jaspar_motif_counts.txt

