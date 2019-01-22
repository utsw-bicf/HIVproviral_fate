#!/bin/bash

#SBATCH --job-name=plotProfile
#SBATCH --partition=super
#SBATCH --output=pP.%j.out
#SBATCH --error=pP.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### This make profile plots
### Do groups for well defined groups, otherwise do individuals

module load deeptools/2.5.0.1

###plot profile
### H3K27ac
bw_loc='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/'
suf="_pooled.fc_signal.bw"
computeMatrix scale-regions \
  -S ${bw_loc}SRR2043614_H3K27ac${suf} ${bw_loc}SRR1603646_H3K27ac${suf} ${bw_loc}SRR1603650_H3K27ac${suf} ${bw_loc}SRR1603654_H3K27ac${suf} ${bw_loc}GSM1603211_GFP_H3K27Ac${suf} \
  -R /project/shared/bicf_workflow_ref/GRCh38/gene.bed \
  --beforeRegionStartLength 3000 \
  --regionBodyLength 5000 \
  --afterRegionStartLength 3000 \
  --skipZeros \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/plotProfile/H3K27ac_matrix.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/plotProfile/H3K27ac_matrix.gz \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/plotProfile/H3K27ac.pdf \
  --perGroup \
  --plotTitle "H3K27ac" \
  --samplesLabel SRR2043614_H3K27ac SRR1603646_H3K27ac SRR1603650_H3K27ac SRR1603654_H3K27ac GSM1603211_GFP_H3K27Ac 

##########
### H3K4m3
bw_loc='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/'
suf="_pooled.fc_signal.bw"
computeMatrix scale-regions \
  -S ${bw_loc}SRR061743_H3K4m3${suf} ${bw_loc}SRR577482_H3K4me3${suf} ${bw_loc}SRR577483_H3K4me3${suf} ${bw_loc}GSM1603213_GFP_H3K4me3${suf} \
  -R /project/shared/bicf_workflow_ref/GRCh38/gene.bed \
  --beforeRegionStartLength 3000 \
  --regionBodyLength 5000 \
  --afterRegionStartLength 3000 \
  --skipZeros \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/plotProfile/H3K4m3_matrix.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/plotProfile/H3K4m3_matrix.gz \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/plotProfile/H3K4m3.pdf \
  --perGroup \
  --plotTitle "H3K4m3" \
  --samplesLabel SRR061743_H3K4m3 SRR577482_H3K4me3 SRR577483_H3K4me3 GSM1603213_GFP_H3K4me3

##########
### each as an individual
samples="SRR2043614_H3K27ac SRR1603646_H3K27ac SRR1603650_H3K27ac SRR1603654_H3K27ac GSM1603211_GFP_H3K27Ac SRR647929_H3K27me3 SRR061743_H3K4m3 SRR577482_H3K4me3 SRR577483_H3K4me3 GSM1603213_GFP_H3K4me3 GSM1603225_GFP_H3K4me1 SRR1522115_med1 SRR3722566_BRD4 SRR3722577_Pol2_DMSO_4h SRR2157609_DMSO_PolII GSM1603217_GFP_PolII GSM1603223_GFP_S5P GSM1603221_GFP_S2P SRR074201_gammaH2AX SRR074195_H2AX SRR074196_H2AX SRR061747_H3K79me2 GSM1603215_GFP_H3K79me3 SRR6010201_YY1 GSM1603209_GFP_H3K36me3 GSM1603227_GFP_H3K9me3"
bw_loc='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/callPeaksMACS/'
suf="_pooled.fc_signal.bw"

for sample in $samples; do
computeMatrix scale-regions \
  -S ${bw_loc}${sample}${suf} \
  -R /project/shared/bicf_workflow_ref/GRCh38/gene.bed \
  --beforeRegionStartLength 3000 \
  --regionBodyLength 5000 \
  --afterRegionStartLength 3000 \
  --skipZeros \
  -o /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/plotProfile/${sample}_matrix.gz

plotProfile -m /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/plotProfile/${sample}_matrix.gz \
  -out /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/plotProfile/${sample}.pdf \
  --perGroup \
  --plotTitle "${sample}" \
  --samplesLabel ${sample}
done
