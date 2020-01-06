#! /bin/bash
########## This script takes the combined output from pastis
########## Draws the Jmol graphic

########## Won't work since there are too many dots; but gave me an idea
########## Plot HIV integration site coordinates vs center dot to see if there is a pattern


########## Add connectors at the end of the file; can only go up to 9999
rm test_connect.txt
for i in {0..9998}
  do 
    j=$((i+1))
    printf "CONECT %4d %4d\n" "${i}" "${j}" >>test_connect.txt
 done

########## Add HIV insertion sites; high, med, low
awk '{OFS = "\t"} $5 <= -0.5 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression_sorted.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_low.bed
awk '{OFS = "\t"} $5 > -0.5 && $5 < 0.5 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression_sorted.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_med.bed
awk '{OFS = "\t"} $5 >= 0.5 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression_sorted.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_high.bed

perl ~/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Chromatin_organization/label_pastis_4_Jmol_HighMedLow.pl \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_high.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_med.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_low.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/combine_20000_abs.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/MDS.HiC_combine_20000.pdb

########## Cat files together
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_added_HighMedLow.pdb /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/test_connect.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_added_HighMedLow_connected.pdb

########## Draw figure
### java -jar /work/BICF/s185797/programs/Jmol/jmol-14.29.55/Jmol.jar
### In console:
# $ color background white;
# $ select :a; color translucent [255, 255, 255];
# $ select ALA; color translucent [230,10,10];spacefill 1
# $ spacefill on;
# $ select ASP; color translucent [20,90,255];spacefill 1
# $ spacefill on;
# $ select SER; color translucent [250,150,0];spacefill 1
# $ spacefill on;
# $ select EDG; color translucent [255, 255, 255];spacefill 0.2
# $ spacefill off;
# $ spacefill on;

# High is red, # Medium is blue, # Low is oragne



########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########## Add HIV insertion sites; high, med, low
awk '{OFS = "\t"} $5 <= -0.5 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression_sorted.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_low.bed
awk '{OFS = "\t"} $5 > -0.5 && $5 < 0.5 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression_sorted.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_med.bed
awk '{OFS = "\t"} $5 >= 0.5 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression_sorted.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_high.bed

perl ~/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Chromatin_organization/label_pastis_4_Jmol_HighMedLow.pl \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_high.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_med.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_low.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/combine_20000_abs.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/MDS.HiC_combine_20000.pdb

head -n 12449 /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_added_HighMedLow.pdb | grep -v "nan" >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_added_HighMedLow_chr1.pdb

grep -v "nan" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_added_HighMedLow.pdb >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_added_HighMedLow_nonan.pdb

########## Draw figure
### java -jar /work/BICF/s185797/programs/Jmol/jmol-14.29.55/Jmol.jar
### In console:
# $ color background white;
# $ select :a; color translucent [255, 255, 255];
# $ select ALA; color translucent [230,10,10];spacefill 1
# $ spacefill on;
# $ select ASP; color translucent [20,90,255];spacefill 1
# $ spacefill on;
# $ select SER; color translucent [250,150,0];spacefill 1
# $ spacefill on;
# $ select EDG; color translucent [255, 255, 255];spacefill 0.2
# $ spacefill off;
# $ spacefill on;

# High is red, # Medium is blue, # Low is oragne


perl ~/Desktop/Holly_git/collaborations/issue141_DOrsoIvan/scripts/Chromatin_organization/label_pastis_4_Jmol.pl \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression_sorted.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/combine_20000_abs.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/MDS.HiC_combine_20000.pdb


########## Grab chr1
ChrList=`awk '{print $1}' combine_20000_abs.bed | uniq`
for i in ${ChrList}; do
grep -w $i combine_20000_abs.bed | tail -n 1
done 
# chr1	248940000	248956422	12448
# chr2	242180000	242193529	24558
# chr3	198280000	198295559	34473
# chr4	190200000	190214555	43984
# chr5	181520000	181538259	53061
# chr6	170800000	170805979	61602
# chr7	159340000	159345973	69570
# chrX	156040000	156040895	77373
# chr8	145120000	145138636	84630
# chr9	138380000	138394717	91550
# chr11	135080000	135086622	98305
# chr10	133780000	133797422	104995
# chr12	133260000	133275309	111659
# chr13	114360000	114364328	117378
# chr14	107040000	107043718	122731
# chr15	101980000	101991189	127831
# chr16	90320000	90338345	132348
# chr17	83240000	83257441	136511
# chr18	80360000	80373285	140530
# chr20	64440000	64444167	143753
# chr19	58600000	58617616	146684
# chrY	57220000	57227415	149546
# chr22	50800000	50818468	152087
# chr21	46700000	46709983	154423

head -n 12449 /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_added.pdb >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000/HIV_added_chr1.pdb


########## Draw figure
### java -jar /work/BICF/s185797/programs/Jmol/jmol-14.29.55/Jmol.jar
### In console:
# $ color background white;
# $ select all; wireframe 0.01;
# $ select :a; color translucent [175, 237, 253];
# $ select ALA; define tel selected; color blue; spacefill 0.02;
