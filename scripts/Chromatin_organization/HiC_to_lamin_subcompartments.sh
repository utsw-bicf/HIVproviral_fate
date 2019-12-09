#! /bin/sh
### Run Hi-C data through HOMER

#SBATCH --job-name=lamin
#SBATCH -p 256GB 
#SBATCH --mem 253952
#SBATCH --output=lamin.%j.out
#SBATCH --error=lamin.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

########## Assign Hi-C results to Lamin A/B; and split A/B compartment

########## Calculate first eigenvector; was already done in standard_Hi-C_HOMER.sh
# ### Compartment analysis
# ### May need to do this again after runnning chromHMM to find open and closed regions
# #echo "Compartment analsysis - PCA"
# export TMPDIR=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/tmp
# runHiCpca.pl auto \
#   /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1 \
#   -cpu 10
# export TMPDIR=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/tmp
# runHiCpca.pl auto \
#   /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2 \
#   -cpu 10
# ###combine outputs
# echo "combine PCA outputs"
# annotatePeaks.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2/lib2.50x100kb.PC1.txt \
#   hg38 \
#   -noblanks \
#   -bedGraph /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2/lib2.50x100kb.PC1.bedGraph /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1/lib1.50x100kb.PC1.bedGraph \
#   > /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/combine_PC1_output.txt
# echo "Compartment analsysis - PCA"

########## Combine tag directories
echo "Combine tag directories"
makeTagDirectory /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged -d /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1 /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2 -tbp 1

########## Calculate PCA analysis
echo "Calculate PCA"
export TMPDIR=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/tmp
runHiCpca.pl auto \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged \
  -cpu 10

########## Find compartments
echo "Find compartments"
findHiCCompartments.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged.PC1.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1_lib2_merged/lib1_lib2_merged_compartments.txt
