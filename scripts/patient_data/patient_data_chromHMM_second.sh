#! /bin/bash

### This makes the pictures for the publication

########## Download E045, E039, E040, and E041
##### For final figures, only use E40
# Website: https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final
# Downloaded to: /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles
# Dowload: E0*_15_coreMarks_hg38lift_segments.bed.gz

# Unzip
zcat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E045_15_coreMarks_hg38lift_segments.bed.gz >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E045_15_coreMarks_hg38lift_segments.bed
zcat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E041_15_coreMarks_hg38lift_segments.bed.gz >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E041_15_coreMarks_hg38lift_segments.bed
zcat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E040_15_coreMarks_hg38lift_segments.bed.gz >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E040_15_coreMarks_hg38lift_segments.bed
zcat /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E039_15_coreMarks_hg38lift_segments.bed.gz >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E039_15_coreMarks_hg38lift_segments.bed

#################################### Create new figures for patient data
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
##### Figure A: Combined analysis; order: BHIVE, Han, Ikeda, Maldarelli, Wagner, Einkauf, Kok, Sharaf, Coffin, McManus, Mack, Garcia-Broncano, Sherril-Mix, Sunshine, Ferris, Lucic
# Put in folder /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient
cp /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/BHIVE.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient/A_BHIVE.bed
cp /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Han2004.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient/B_Han2004.bed
cp /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Ikeda2007.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient/C_Ikeda2007.bed
cp /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Maldarelli2014.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient/D_Maldarelli2014.bed
cp /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Wagner2014.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient/E_Wagner2014.bed
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Einkauf2019.bed  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Einkauf_patient1.bed  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Einkauf_patient2.bed  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Einkauf_patient3.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient/F_Einkauf.bed
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Kok_2016/Kok_2016_patient*_CD4_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient/G_Kok_2016.bed
cp /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Sharaf2018.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient/H_Sharaf2018.bed
cp /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Coffin2019.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient/I_Coffin2019.bed
cp /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/McManus2019.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient/J_McManus2019.bed
cp /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Mack2003.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient/K_Mack2003.bed
cp /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Garcia-Broncano2019.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient/L_Garcia-Broncano2019.bed
cp /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Sherrill-Mix2013.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient/M_Sherrill-Mix2013.bed
cp /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Sunshine2016.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient/N_Sunshine2016.bed
cp /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Ferris2019.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient/O_Ferris2019.bed
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Lusic2019/genome_structure_and_HIV_integration/GSE134382_NIH/Lusic2019_*_hg38_liftover.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient/P_Lusic_2019.bed

### Do overlap plot for chromHMM
## Jurkat
#unset DISPLAY
#java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
#  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/Jurkat_CombinedPatient

##E045
#unset DISPLAY
#java -Djava.awt.headless=true -Xmx1000g -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
#  -b 1 \
#  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/labelmappingfile_known.txt \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E045_15_coreMarks_hg38lift_segments.bed \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/E045_CombinedPatient

##E039
#unset DISPLAY
#java -Djava.awt.headless=true -Xmx1000g -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
#  -b 1 \
#  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/labelmappingfile_known.txt \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E039_15_coreMarks_hg38lift_segments.bed \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/E039_CombinedPatient

#E040
unset DISPLAY
java -Djava.awt.headless=true -Xmx1000g -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -b 1 \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/labelmappingfile_known.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E040_15_coreMarks_hg38lift_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/E040_CombinedPatient

##E041
#unset DISPLAY
#java -Djava.awt.headless=true -Xmx1000g -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
#  -b 1 \
#  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/labelmappingfile_known.txt \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E041_15_coreMarks_hg38lift_segments.bed \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_CombinedPatient \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/E041_CombinedPatient


########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
##### Figure B: Individual patient analysis: Kok and Einkauf
# Put in folder /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IndividualPatient

for x in /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/Kok_2016/Kok_2016_patient*_CD4_hg38_liftover.bed; do
  header=$(basename ${x} _CD4_hg38_liftover.bed)
  cp ${x} /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IndividualPatient/${header}.bed
done

cp /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Einkauf_patient1.bed  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Einkauf_patient2.bed  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Einkauf_patient3.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IndividualPatient


### Do overlap plot for chromHMM
# Jurkat
unset DISPLAY
java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IndividualPatient \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureB_IndividualPatient/Jurkat_IndividualPatient

#E045
unset DISPLAY
java -Djava.awt.headless=true -Xmx1000g -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -b 1 \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/labelmappingfile_known.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E045_15_coreMarks_hg38lift_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IndividualPatient\
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureB_IndividualPatient/E045_IndividualPatient

#E039
unset DISPLAY
java -Djava.awt.headless=true -Xmx1000g -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -b 1 \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/labelmappingfile_known.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E039_15_coreMarks_hg38lift_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IndividualPatient \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureB_IndividualPatient/E039_IndividualPatient

#E040
unset DISPLAY
java -Djava.awt.headless=true -Xmx1000g -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -b 1 \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/labelmappingfile_known.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E040_15_coreMarks_hg38lift_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IndividualPatient \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureB_IndividualPatient/E040_IndividualPatient

#E041
unset DISPLAY
java -Djava.awt.headless=true -Xmx1000g -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -b 1 \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/labelmappingfile_known.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E041_15_coreMarks_hg38lift_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IndividualPatient \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureB_IndividualPatient/E041_IndividualPatient




########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
##### Figure C: Intact vs Defective; Einkauf
# Put in folder /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IntactDefective

grep "Intact" /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Einkauf_patient1.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IntactDefective/A_Einkauf_patient1_Intact.bed
grep -v "Intact" /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Einkauf_patient1.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IntactDefective/D_Einkauf_patient1_Defective.bed

grep "Intact" /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Einkauf_patient2.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IntactDefective/B_Einkauf_patient2_Intact.bed
grep -v "Intact" /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Einkauf_patient2.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IntactDefective/E_Einkauf_patient2_Defective.bed

grep "Intact" /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Einkauf_patient3.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IntactDefective/C_Einkauf_patient3_Intact.bed
grep -v "Intact" /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Einkauf_patient3.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IntactDefective/F_Einkauf_patient3_Defective.bed

grep "Intact" /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Garcia-Broncano2019.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IntactDefective/G_Garcia-Broncano_Intact.bed
grep -v "Intact" /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_input/Garcia-Broncano2019.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IntactDefective/H_Garcia-Broncano_Defective.bed

### Do overlap plot for chromHMM
## Jurkat
#unset DISPLAY
#java -Djava.awt.headless=true -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
#  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/labelmappingfile.txt \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/jurkat_15_segments_reorder.bed \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IntactDefective \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureC_IntactDefective/Jurkat_IntactDefective

##E045
#unset DISPLAY
#java -Djava.awt.headless=true -Xmx1000g -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
#  -b 1 \
#  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/labelmappingfile_known.txt \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E045_15_coreMarks_hg38lift_segments.bed \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IntactDefective\
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureC_IntactDefective/E045_IntactDefective

##E039
#unset DISPLAY
#java -Djava.awt.headless=true -Xmx1000g -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
#  -b 1 \
#  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/labelmappingfile_known.txt \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E039_15_coreMarks_hg38lift_segments.bed \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IntactDefective \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureC_IntactDefective/E039_IntactDefective

#E040
unset DISPLAY
java -Djava.awt.headless=true -Xmx1000g -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
  -b 1 \
  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/labelmappingfile_known.txt \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E040_15_coreMarks_hg38lift_segments.bed \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IntactDefective \
  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureC_IntactDefective/E040_IntactDefective

##E041
#unset DISPLAY
#java -Djava.awt.headless=true -Xmx1000g -jar /work/BICF/s185797/programs/ChromHMM/ChromHMM.jar OverlapEnrichment \
#  -b 1 \
#  -m /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureA_CombinedPatient/labelmappingfile_known.txt \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/known_state_bedfiles/E041_15_coreMarks_hg38lift_segments.bed \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/Patient_data/overlap_IntactDefective \
#  /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/output_histone7/histone_learn_15/Reorder_sameasErnst/FigureC_IntactDefective/E041_IntactDefective
