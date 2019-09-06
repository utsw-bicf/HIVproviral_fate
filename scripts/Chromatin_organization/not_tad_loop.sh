#! /bin/bash

########## This script creates a bed file that isn't a tad or isn't a loop
module load bedops/2.4.14

### Fix loop file coordinates
awk '{OFS="\t"}; $1==$4 { print $1, $2, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/Loop_outercoord.bed

bedops --complement /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/Loop_outercoord.bed > /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/not_loop.bed

### Manually add each chr start and end to file


##########
### Fix tad file coordinates
awk '{OFS="\t"}; $1==$4 { print $1, $2, $3}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.tad.2D.bed | sort -k 1,1 -k 2,2n >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/Tad_coord.bed

bedops --complement /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/Tad_coord.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/not_tad.bed

### Manually add each chr start and end to file
