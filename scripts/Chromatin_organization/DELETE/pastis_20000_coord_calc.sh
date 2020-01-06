#! /bin/bash
########## This script was an idea after plotting didn't work
########## Take HIV integration High, Medium, Low and
########## Take the x,y,z coordinates and calculate distance to center

wd='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HiCpro/output_111419/hic_results/matrix/combine/raw/20000'

##### 1. Find center
### Min and Max of x; x = 0.01, 199.997
grep -v "nan" ${wd}/MDS.HiC_combine_20000.pdb | sort -k 6,6n | head -n 1
grep -v "nan" ${wd}/MDS.HiC_combine_20000.pdb | sort -k 6,6n | tail -n 1

### Min and Max of y; y = 0, 199.999
grep -v "nan" ${wd}/MDS.HiC_combine_20000.pdb | sort -k 7,7n | head -n 1
grep -v "nan" ${wd}/MDS.HiC_combine_20000.pdb | sort -k 7,7n | tail -n 1 

### Min and Max of z; z = 0.02, 199.999
grep -v "nan" ${wd}/MDS.HiC_combine_20000.pdb | sort -k 8,8n | head -n 1
grep -v "nan" ${wd}/MDS.HiC_combine_20000.pdb | sort -k 8,8n | tail -n 1 

# So center is 100,100,100


##### 2. Separate HIV into low, medium, and high with cutoffs of -0.5 and 0.5
awk '{OFS = "\t"} $5 <= -0.5 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression_sorted.bed >${wd}/HIV_low.bed
awk '{OFS = "\t"} $5 > -0.5 && $5 < 0.5 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression_sorted.bed >${wd}/HIV_med.bed
awk '{OFS = "\t"} $5 >= 0.5 {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression_sorted.bed | grep -v "chrUn" >${wd}/HIV_high.bed

# High = 633, Medium = 794, Low = 134


##### 3. Get bin of each HIV
module load bedtools
bedtools intersect -wao -a ${wd}/HIV_low.bed -b ${wd}/combine_20000_abs.bed >${wd}/HIV_low_bin.bed
bedtools intersect -wao -a ${wd}/HIV_med.bed -b ${wd}/combine_20000_abs.bed >${wd}/HIV_med_bin.bed
bedtools intersect -wao -a ${wd}/HIV_high.bed -b ${wd}/combine_20000_abs.bed >${wd}/HIV_high_bin.bed


##### 4. Match HIV bin with 20,000 bin x,y,z coordinates
##### 5. Add distance to center sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1) ^2))
module load R
Rscript merge_HIVbins_pdb.R


##### 6. Find min, max, avg, mean of distances; done in R
# > min(pdb_high$dist)
# [1] 15.30979
# > max(pdb_high$dist)
# [1] 168.4169
# > mean(pdb_high$dist)
# [1] 98.23211
# > median(pdb_high$dist)
# [1] 100.8168

# > min(pdb_med$dist)
# [1] 14.33311
# > max(pdb_med$dist)
# [1] 159.3368
# > mean(pdb_med$dist)
# [1] 96.38333
# > median(pdb_med$dist)
# [1] 97.56887


# > min(pdb_low$dist)
# [1] 24.35604
# > max(pdb_low$dist)
# [1] 165.7815
# > mean(pdb_low$dist)
# [1] 97.9439
# > median(pdb_low$dist)
# [1] 102.229
