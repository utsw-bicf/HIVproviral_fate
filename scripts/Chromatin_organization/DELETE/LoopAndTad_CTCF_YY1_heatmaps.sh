#! /bin/bash
#SBATCH --job-name=hm
#SBATCH -p 256GB
#SBATCH --mem 252928
#SBATCH --output=hm.%j.out
#SBATCH --error=hm.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END


########## This script does a heatmap of CTCF and YY1 on Tad and Loop boundaries

########## CTCF
CTCF='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/chipseq_analysis_4CTCF/workflow/output_CTCF/callPeaksMACS/CTCF_YOUNG_GSM1689152_pooled.fc_signal.bw'

########## YY1
##### Merge bam files, and remove dups
#module load samtools/1.6
#samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/YY1_merged.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/output_031419/bowtie_results/bwt2/rep1/SRR6010260_GRCh38.bwt2pairs.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/output_031419/bowtie_results/bwt2/rep2/SRR6010261_GRCh38.bwt2pairs.bam
#samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/YY1_merged.bam

#module load sambamba/0.6.6
#sambamba markdup -r -t 32 --hash-table-size=17592186044416 --overflow-list-size=20000000 --io-buffer-size=256 /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/YY1_merged.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/YY1_merged_dups.bam
#samtools sort -@ 32 /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/YY1_merged_dups.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/YY1.bam
#samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/YY1.bam


##### Make bw file
#module load deeptools/2.5.0.1
#bamCoverage --numberOfProcessors=max/2 -b /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/YY1.bam -o /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/YY1.bw

YY1='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/YY1.bw'


########## Loop
##### Create a bed file of left and right loop boundaries
#grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D.bed | awk '{OFS = "\t"} {print $1, $2, $3}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/left_loop.bed
ll='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/left_loop.bed'

#grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D.bed | awk '{OFS = "\t"} {print $4, $5, $6}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/right_loop.bed
rl='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/right_loop.bed'


########## Make heatmap
module load deeptools/2.5.0.1
computeMatrix scale-regions -S ${CTCF} \
  -R ${ll} ${rl} \
  -a 15000 \
  -b 15000 \
  -m 3500 \
  -p max/2 \
  -o CTCF_loop.mat.gz

plotHeatmap -m CTCF_loop.mat.gz --startLabel LB --endLabel RB --samplesLabel CTCF --regionsLabel Left_Loop Right_loop --colorList 'white,black' -out CTCF_loop.pdf


module load deeptools/2.5.0.1
computeMatrix scale-regions -S ${YY1} \
  -R ${ll} ${rl} \
  -a 15000 \
  -b 15000 \
  -m 3500 \
  -p max/2 \
  -o YY1_loop.mat.gz

plotHeatmap -m YY1_loop.mat.gz --startLabel LB --endLabel RB --samplesLabel YY1 --regionsLabel Left_Loop Right_loop --colorList 'white,black' -out YY1_loop.pdf


########## Tad
grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.tad.2D.bed | awk '{OFS = "\t"} {print $1, $2, $2+1}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/left_tad.bed
lt='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/left_tad.bed'
grep -v "#" /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.tad.2D.bed | awk '{OFS = "\t"} {print $1, $3-1, $3}' >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/right_tad.bed
rt='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/heatmaps_boundaries/right_tad.bed'

computeMatrix reference-point -S ${CTCF} \
  -R ${lt} ${rt} \
  -a 15000 \
  -b 15000 \
  --referencePoint center \
  -p max/2 \
  -o CTCF_tad.mat.gz

plotHeatmap -m CTCF_tad.mat.gz --refPointLabel "Tad boundary" --samplesLabel CTCF --regionsLabel Left_Tad Right_Tad --colorList 'white,blue' -out CTCF_tad.pdf


computeMatrix reference-point -S ${YY1} \
  -R ${lt} ${rt} \
  -a 15000 \
  -b 15000 \
  --referencePoint center \
  -p max/2 \
  -o YY1_tad.mat.gz

plotHeatmap -m YY1_tad.mat.gz --refPointLabel "Tad boundary" --samplesLabel YY1 --regionsLabel Left_Tad Right_Tad --colorList 'white,blue' -out YY1_tad.pdf


########## Do H3K27ac just for reference
H3K27ac='/project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/callPeaksMACS/H3K27ac_pooled.fc_signal.bw'


computeMatrix scale-regions -S ${H3K27ac} \
  -R ${ll} ${rl} \
  -a 15000 \
  -b 15000 \
  -m 3500 \
  -p max/2 \
  -o H3K27ac_loop.mat.gz

plotHeatmap -m H3K27ac_loop.mat.gz --startLabel LB --endLabel RB --samplesLabel H3K27ac --regionsLabel Left_Loop Right_loop --colorList 'white,black' -out H3K27ac_loop.pdf


computeMatrix reference-point -S ${H3K27ac} \
  -R ${lt} ${rt} \
  -a 15000 \
  -b 15000 \
  --referencePoint center \
  -p max/2 \
  -o H3K27ac_tad.mat.gz

plotHeatmap -m H3K27ac_tad.mat.gz --refPointLabel "Tad boundary" --samplesLabel H3K27ac --regionsLabel Left_Tad Right_Tad --colorList 'white,blue' -out H3K27ac_tad.pdf

