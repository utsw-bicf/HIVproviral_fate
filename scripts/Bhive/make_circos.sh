#! /bin/bash
### This script makes the circos plots for each group 1-6
### run inside of folder: /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/circos
### Each group needs their own group*.conf for colors
### group 1 = blue
### group 2 = red
### group 3 = orange
### group 4 = green
### group 5 = grey
### group 6 = purple

module load circos/0.696
module load bedtools/2.25.0

group1="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group1_Intergenic_same.bed"
group2="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group2_Intergenic_downstream_opposite.bed"
group3="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group3_Intergenic_upstream_opposite.bed"
group4="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group4_Intragenic_same_1gene.bed"
group5="/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group5_Intragenic_opposite_1gene.bed"

### Fix group 6,7,8
bedtools intersect -wao -a <(awk '{OFS="\t"} {print $1,$2-2, $3+2}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6_2gene_opposite_noexp.bed) \
  -b /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.bed \
  >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6.bed
bedtools intersect -wao -a <(awk '{OFS="\t"} {print $1,$2-2, $3+2}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group7_2gene_same_noexp.bed) \
  -b /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.bed \
  >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group7.bed
bedtools intersect -wao -a <(awk '{OFS="\t"} {print $1,$2-2, $3+2}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group8_2gene_opposite_noexp.bed) \
  -b /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.bed \
  >/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group8.bed



########## Split each group by positive and negitive;change chr to hs
# Group 1
awk '{OFS="\t"} $11 >= 0 {print $1, $2, $3, $11}' ${group1} | perl -p -e 's|^chr|hs|g' >g1pos.txt
awk '{OFS="\t"} $11 < 0 {print $1, $2, $3, $11}' ${group1} | perl -p -e 's|^chr|hs|g' >g1neg.txt

circos -conf group1_circos.conf -outputfile group1 -debug_group summary,timer >run_g1.out

# Group 2
awk '{OFS="\t"} $11 >= 0 {print $1, $2, $3, $11}' ${group2} | perl -p -e 's|^chr|hs|g' >g2pos.txt
awk '{OFS="\t"} $11 < 0 {print $1, $2, $3, $11}' ${group2} | perl -p -e 's|^chr|hs|g' >g2neg.txt

circos -conf group2_circos.conf -outputfile group2 -debug_group summary,timer >run_g2.out

# Group 3
awk '{OFS="\t"} $11 >= 0 {print $1, $2, $3, $11}' ${group3} | perl -p -e 's|^chr|hs|g' >g3pos.txt
awk '{OFS="\t"} $11 < 0 {print $1, $2, $3, $11}' ${group3} | perl -p -e 's|^chr|hs|g' >g3neg.txt

circos -conf group3_circos.conf -outputfile group3 -debug_group summary,timer >run_g3.out

# Group 4
awk '{OFS="\t"} $11 >= 0 {print $1, $2, $3, $11}' ${group4} | perl -p -e 's|^chr|hs|g' >g4pos.txt
awk '{OFS="\t"} $11 < 0 {print $1, $2, $3, $11}' ${group4} | perl -p -e 's|^chr|hs|g' >g4neg.txt

circos -conf group4_circos.conf -outputfile group4 -debug_group summary,timer >run_g4.out


# Group 5
awk '{OFS="\t"} $11 >= 0 {print $1, $2, $3, $11}' ${group5} | perl -p -e 's|^chr|hs|g' >g5pos.txt
awk '{OFS="\t"} $11 < 0 {print $1, $2, $3, $11}' ${group5} | perl -p -e 's|^chr|hs|g' >g5neg.txt

circos -conf group5_circos.conf -outputfile group5 -debug_group summary,timer >run_g5.out


# Group 6
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group7.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group8.bed | awk '{OFS="\t"} $8 >= 0 {print $4, $5, $6, $8}' | perl -p -e 's|^chr|hs|g' >g6pos.txt
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group6.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group7.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/expression_annotated_bedfiles/group8.bed | awk '{OFS="\t"} $8 < 0 {print $4, $5, $6, $8}' | perl -p -e 's|^chr|hs|g' >g6neg.txt

circos -conf group6_circos.conf -outputfile group6 -debug_group summary,timer >run_g6.out
