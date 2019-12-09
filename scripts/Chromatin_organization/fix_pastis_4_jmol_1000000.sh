#! /bin/bash
########## This script fixes pastis output for Jmol

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
rm test_connect.txt
for i in {0..3101}
  do 
    j=$((i+1))
    printf "CONECT %4d %4d\n" "${i}" "${j}" >>test_connect.txt
 done

cat MDS.HiC_combine_1000000.pdb test_connect.txt >test4.pdb

grep -v "CONECT  249  250" test4.pdb | grep -v "CONECT  492  493" | grep -v "CONECT  691  692" | grep -v "CONECT  882  883" | grep -v "CONECT 1064 1065" | grep -v "CONECT 1235 1236" | grep -v "CONECT  1395  1396" | grep -v "CONECT  1552  1553" | grep -v "CONECT  1698  1699" |grep -v "CONECT 1837 1838" |grep -v "CONECT 1973 1974" |grep -v "CONECT 2107 2108" |grep -v "CONECT 2241 2242" |grep -v "CONECT 2356 2357" |grep -v "CONECT 2464 2465" |grep -v "CONECT 2566 2567" |grep -v "CONECT 2657 2658" |grep -v "CONECT 2741 2742" |grep -v "CONECT 2822 2823" |grep -v "CONECT 2887 2888" |grep -v "CONECT 2946 2947" | grep -v "CONECT 3004 3005" | grep -v "CONECT 3055 3056" | grep -v "CONECT 3102 3103" >test5.pdb

### Relabel chromosomes; must keep same tabbing;
# Do by head/tail for each chromosome, change column 5, cat chromosomes back together

### add HIV
perl ~/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Chromatin_organization/label_pastis_4_Jmol.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression_sorted.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/1000000/combine_1000000_abs.bed test5.pdb

### java -jar /work/BICF/s185797/programs/Jmol/jmol-14.29.55/Jmol.jar
### In console:
# $ color background white;
# $ select all; wireframe 0.01;
# $ select :a; color translucent [175, 237, 253];
# $ select ALA; define tel selected; color blue; spacefill 0.02;

########## Below is the info used to determine chromosomes
# ChrList=`awk '{print $1}' combine_1000000_abs.bed | uniq`
# for i in ${ChrList}; do
#  grep -w $i combine_1000000_abs.bed | tail -n 1
# done 

#chr1	248000000	248956422	249
#chr2	242000000	242193529	492
#chr3	198000000	198295559	691
#chr4	190000000	190214555	882
#chr5	181000000	181538259	1064
#chr6	170000000	170805979	1235
#chr7	159000000	159345973	1395
#chrX	156000000	156040895	1552
#chr8	145000000	145138636	1698
#chr9	138000000	138394717	1837
#chr11	135000000	135086622	1973
#chr10	133000000	133797422	2107
#chr12	133000000	133275309	2241
#chr13	114000000	114364328	2356
#chr14	107000000	107043718	2464
#chr15	101000000	101991189	2566
#chr16	90000000	90338345	2657
#chr17	83000000	83257441	2741
#chr18	80000000	80373285	2822
#chr20	64000000	64444167	2887
#chr19	58000000	58617616	2946
#chrY	57000000	57227415	3004
#chr22	50000000	50818468	3055
#chr21	46000000	46709983	3102


