#! /bin/bash

########## This script tells where to download patient data,
########## Organize patient data,
########## and make chromHMM

########## Table: /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Table_Figure8_patientdata.docx


########## 1) RID: retrovirus integration database; 
##### Downloaded from https://rid.ncifcrf.gov/; but is in hg19
### Lift over, make bed file and use https://genome.ucsc.edu/cgi-bin/hgLiftOver
awk '{OFS = "\t"} NR > 1 {print $2, $3-1, $3, $1, $6, $5, $NF}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/1130815007.dwnld.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg19.bed

module load UCSC_userApps/v317
liftOver /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg19.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/hg19ToHg38.over.chain.gz \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_liftover.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_notlifted.bed

### Separate by pubmedID
awk '$7 ~ /12843741/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_pubmed12843741.bed
awk '$7 ~ /15163705/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_pubmed15163705.bed
awk '$7 ~ /17262715/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_pubmed17262715.bed
awk '$7 ~ /23953889/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_pubmed23953889.bed
awk '$7 ~ /24968937/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_pubmed24968937.bed
awk '$7 ~ /25011556/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_pubmed25011556.bed
awk '$7 ~ /26545813/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_pubmed26545813.bed
awk '$7 ~ /26912621/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_pubmed26912621.bed
awk '$7 ~ /30024859/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_pubmed30024859.bed
awk '$7 ~ /30688658/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_pubmed30688658.bed
awk '$7 ~ /30857886/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_pubmed30857886.bed
awk '$7 ~ /31217357/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_pubmed31217357.bed
awk '$7 ~ /31291371/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_pubmed31291371.bed
awk '$7 ~ /31361603/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_pubmed31361603.bed
awk '$7 ~ /SRP065157/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/RID_retrovirus_integration_database/RID_hg38_pubmedSRP065157.bed


########## 2) Einkauf 2019; Got directly from supplementary table
##### Make excel into txt file for each patient, then make bed
awk '{OFS = "\t"} NR > 1 {print $3, $4-1, $4, "patient1", $2, $5}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Einkauf_2019/Einkauf_patient1.txt | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Einkauf_2019/Einkauf_patient1.bed
awk '{OFS = "\t"} NR > 1 {print $3, $4-1, $4, "patient1", $2, $5}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Einkauf_2019/Einkauf_patient2.txt | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Einkauf_2019/Einkauf_patient2.bed
awk '{OFS = "\t"} NR > 1 {print $3, $4-1, $4, "patient1", $2, $5}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Einkauf_2019/Einkauf_patient3.txt | sort -k 1,1 -k 2,2n | uniq >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Einkauf_2019/Einkauf_patient3.bed


########## 3) Brady 2009**, (Bushman lab); Not found online

########## 4) Wagner 2014**, (Frenkel lab); RID pubmed ID: 25011556

########## 5) Maldarelli 2014**, (Hughes lab); RID pubmed ID: 24968937

########## 6) Han 2004**, (Siliciano lab); RID pubmed ID: 15163705

########## 7) Cohn 2015**, Nussenzweig lab; Not found online

########## 8) Ikeda 2007**, (Matsushita lab); RID pubmed ID: 15163705

########## 9) Kok 2016**,(Metzner lab); From supplementary table S4
##### Convert pdf to txt with Tabula
## Make bed file
awk '{OFS = "\t"} NR >1 && $1 == 1 {print $6, $4-1, $4, $5, $2, $7}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Kok_2016/Kok_2016.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Kok_2016/Kok_2016_patient1.bed
awk '{OFS = "\t"} NR >1 && $1 == 2 {print $6, $4-1, $4, $5, $2, $7}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Kok_2016/Kok_2016.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Kok_2016/Kok_2016_patient2.bed
awk '{OFS = "\t"} NR >1 && $1 == 3 {print $6, $4-1, $4, $5, $2, $7}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Kok_2016/Kok_2016.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Kok_2016/Kok_2016_patient3.bed
awk '{OFS = "\t"} NR >1 && $1 == 4 {print $6, $4-1, $4, $5, $2, $7}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Kok_2016/Kok_2016.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Kok_2016/Kok_2016_patient4.bed
awk '{OFS = "\t"} NR >1 && $1 == 5 {print $6, $4-1, $4, $5, $2, $7}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Kok_2016/Kok_2016.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Kok_2016/Kok_2016_patient5.bed
awk '{OFS = "\t"} NR >1 && $1 == 6 {print $6, $4-1, $4, $5, $2, $7}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Kok_2016/Kok_2016.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Kok_2016/Kok_2016_patient6.bed
awk '{OFS = "\t"} NR >1 && $1 == 7 {print $6, $4-1, $4, $5, $2, $7}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Kok_2016/Kok_2016.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Kok_2016/Kok_2016_patient7.bed


########## 10) Lucic 2019; Got directly from NCBI GSE134382
### Make as bed file
awk -F "," '{OFS = "\t"} NR >1 {print "chr"$1,$3-1,$3,$4, $5, $2}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/GSE134382_33INF_processed_IS.csv | perl -p -e 's|minus|-|g;s|plus|+|g;' >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/Lusic2019_33INF.bed
awk -F "," '{OFS = "\t"} NR >1 {print "chr"$1,$3-1,$3,$4, $5, $2}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/GSE134382_33mockINF_processed_IS.csv | perl -p -e 's|minus|-|g;s|plus|+|g;' >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/Lusic2019_33mockINF.bed
awk -F "," '{OFS = "\t"} NR >1 {print "chr"$1,$3-1,$3,$4, $5, $2}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/GSE134382_34INF_processed_IS.csv | perl -p -e 's|minus|-|g;s|plus|+|g;' >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/Lusic2019_34INF.bed
awk -F "," '{OFS = "\t"} NR >1 {print "chr"$1,$3-1,$3,$4, $5, $2}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/GSE134382_34mockINF_processed_IS.csv | perl -p -e 's|minus|-|g;s|plus|+|g;' >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/Lusic2019_34mockINF.bed
awk -F "," '{OFS = "\t"} NR >1 {print "chr"$1,$3-1,$3,$4, $5, $2}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/GSE134382_3673A_processed_IS.csv | perl -p -e 's|minus|-|g;s|plus|+|g;' >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/Lusic2019_3673A.bed
awk -F "," '{OFS = "\t"} NR >1 {print "chr"$1,$3-1,$3,$4, $5, $2}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/GSE134382_3673Q_processed_IS.csv | perl -p -e 's|minus|-|g;s|plus|+|g;' >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/Lusic2019_3673Q.bed
#awk -F "," '{OFS = "\t"} NR >1 {print "chr"$1,$3-1,$3,$4, $5, $2}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/GSE134382_6612A_processed_IS.csv| perl -p -e 's|minus|-|g;s|plus|+|g;' >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/Lusic2019_6612A.bed
awk '{OFS = "\t"} NR >1 {print "chr"$1,$3-1,$3,$4, $5, $2}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/GSE134382_6612A_processed_IS_fixed.txt | perl -p -e 's|minus|-|g;s|plus|+|g;' >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/Lusic2019_6612A.bed
awk -F "," '{OFS = "\t"} NR >1 {print "chr"$1,$3-1,$3,$4, $5, $2}' /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/GSE134382_6612Q_processed_IS.csv | perl -p -e 's|minus|-|g;s|plus|+|g;' >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/Lusic2019_6612Q.bed

########## 11) HIV Bhive; use
# /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.bed


######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
### Make table in excel of all the data

### Do overlap plot for chromHMM
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/overlap_patient
