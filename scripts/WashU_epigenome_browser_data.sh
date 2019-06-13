#! /bin/sh
### This bash script makes the files for WashU epigenome browser
module load UCSC_userApps/v317
module load bedtools
module load samtools
module load deeptools/2.5.0.1
module load juicebox/1.5.6

cd /archive/shared/DOrso_BICF/WashU_browser_data

### DNAse-seq; replicated narrowPeak to bigBed
awk '{OFS="\t"}; $1 !~ /_/ {print $0}' '/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_DNase/workflow/output/consensusPeaks/DNAse-seq_ENCODE_GSM736501.replicated.narrowPeak' | sort -k 1,1 -k 2,2n >delete.narrowPeak
bedToBigBed -type=bed6+4 delete.narrowPeak /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes DNase-seq_replicated_peaks.bigBed
rm delete.narrowPeak


### MNase-seq; iNPS significant peaks bed format to bigBed
grep -v "1\:Chrom" /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/iNPS/iNPS_output/MNase_significant_all.bed | sort -k 1,1 -k 2,2n | uniq >delete.bed
bedToBigBed -type=bed5+5 delete.bed /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes MNase-seq_iNPS_significant.bigBed
rm delete.bed


### ChiP-seq H3K27ac; replicated narrowPeak from consensus peaks all experiments
awk '{OFS="\t"}; $1 !~ /_/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/consensusPeaks/H3K27ac.replicated.narrowPeak | sort -k 1,1 -k 2,2n >delete.narrowPeak
bedToBigBed -type=bed6+4 delete.narrowPeak /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes H3K27ac_replicated_peaks.bigBed
rm delete.narrowPeak


### ChiP-seq H3K4me3; replicated narrowPeak from consensus peaks all experiment
awk '{OFS="\t"}; $1 !~ /_/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/consensusPeaks/H3K4me3.replicated.narrowPeak | sort -k 1,1 -k 2,2n >delete.narrowPeak
bedToBigBed -type=bed6+4 delete.narrowPeak /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes H3K4me3_replicated_peaks.bigBed
rm delete.narrowPeak


### ChiP-seq H3K79me2; replicated narrowPeak; single experiment
awk '{OFS="\t"}; $1 !~ /_/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/SRR061747_H3K79me2.replicated.narrowPeak | sort -k 1,1 -k 2,2n >delete.narrowPeak
bedToBigBed -type=bed6+4 delete.narrowPeak /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes H3K79me2_replicated_peaks.bigBed
rm delete.narrowPeak


### ChiP-seq H3K27me3; replicated narrowPeak from consensus peaks all experiments
awk '{OFS="\t"}; $1 !~ /_/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/SRR647929_H3K27me3.replicated.narrowPeak | sort -k 1,1 -k 2,2n >delete.narrowPeak
bedToBigBed -type=bed6+4 delete.narrowPeak /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes H3K27me3_replicated_peaks.bigBed
rm delete.narrowPeak


### ChiP-seq PolII; replicated narrowPeak from consensus peaks all experiment
awk '{OFS="\t"}; $1 !~ /_/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus_PolII/consensusPeaks/PolII.replicated.narrowPeak | sort -k 1,1 -k 2,2n >delete.narrowPeak
bedToBigBed -type=bed6+4 delete.narrowPeak /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes PolII_replicated_peaks.bigBed
rm delete.narrowPeak


### ChiP-seq H3K4me1; replicated narrowPeak; single experiment
awk '{OFS="\t"}; $1 !~ /_/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/GSM1603225_GFP_H3K4me1.replicated.narrowPeak | sort -k 1,1 -k 2,2n >delete.narrowPeak
bedToBigBed -type=bed6+4 delete.narrowPeak /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes H3K4me1_replicated_peaks.bigBed
rm delete.narrowPeak


### ChiP-seq H3K36me3; replicated narrowPeak; single experiment
awk '{OFS="\t"}; $1 !~ /_/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/GSM1603209_GFP_H3K36me3.replicated.narrowPeak | sort -k 1,1 -k 2,2n >delete.narrowPeak
bedToBigBed -type=bed6+4 delete.narrowPeak /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes H3K36me3_replicated_peaks.bigBed
rm delete.narrowPeak


### ChiP-seq H3K79me3; replicated narrowPeak; single experiment
awk '{OFS="\t"}; $1 !~ /_/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/GSM1603215_GFP_H3K79me3.replicated.narrowPeak | sort -k 1,1 -k 2,2n >delete.narrowPeak
bedToBigBed -type=bed6+4 delete.narrowPeak /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes H3K79me3_replicated_peaks.bigBed
rm delete.narrowPeak


### ChiP-seq H3K9me3; replicated narrowPeak; single experiment
awk '{OFS="\t"}; $1 !~ /_/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/GSM1603227_GFP_H3K9me3.replicated.narrowPeak | sort -k 1,1 -k 2,2n >delete.narrowPeak
bedToBigBed -type=bed6+4 delete.narrowPeak /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes H3K9me3_replicated_peaks.bigBed
rm delete.narrowPeak


### ChiP-seq BRD4; replicated narrowPeak; single experiment
awk '{OFS="\t"}; $1 !~ /_/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/SRR3722566_BRD4.replicated.narrowPeak | sort -k 1,1 -k 2,2n >delete.narrowPeak
bedToBigBed -type=bed6+4 delete.narrowPeak /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes BRD4_replicated_peaks.bigBed
rm delete.narrowPeak


### ChiP-seq med1; replicated narrowPeak; single experiment
awk '{OFS="\t"}; $1 !~ /_/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/SRR1522115_med1.replicated.narrowPeak | sort -k 1,1 -k 2,2n >delete.narrowPeak
bedToBigBed -type=bed6+4 delete.narrowPeak /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes med1_replicated_peaks.bigBed
rm delete.narrowPeak


### ChiP-seq YY1; replicated narrowPeak; single experiment
awk '{OFS="\t"}; $1 !~ /_/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/SRR6010201_YY1.replicated.narrowPeak | sort -k 1,1 -k 2,2n >delete.narrowPeak
bedToBigBed -type=bed6+4 delete.narrowPeak /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes YY1_replicated_peaks.bigBed
rm delete.narrowPeak


### tt-seq; bam to bed of merged bam to bigBed
## Done previously
#samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Ttseq_CRAMER_merged.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260187.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_TT/Ttseq_CRAMER_GSM2260188.dedup.bam
bamCoverage -b /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/Ttseq_CRAMER_merged.bam -o ttseq_merged.bw


### RNA-seq; bam to bed of Emily replicates to bigBed
## Done previously
#samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/RNAseq.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq_reg/workflow/output_PE/Emily_jurkat_*.dedup.bam
bamCoverage -b /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/annotate/1D_expression_vs_chip/RNAseq.bam -o RNAseq_Emily_merged.bw


### 4Su-seq; bam to bed of merged bam to bigBed
## Done previously
samtools merge /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_4sU/4sUseq_merged.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_4sU/S4sUseq_KEENE_GSM1833449.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_4sU/S4sUseq_KEENE_GSM1833448.dedup.bam /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_4sU/S4sUseq_KEENE_GSM1833447.dedup.bam
samtools index /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_4sU/4sUseq_merged.bam
bamCoverage -b /project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/workflow/output_4sU/4sUseq_merged.bam -o 4sUseq_merged.bw


### CTCF-seq; replicated narrowPeak; single experiment
awk '{OFS="\t"}; $1 !~ /_/ {print $0}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/chipseq_analysis_4CTCF/workflow/output_CTCF/consensusPeaks/CTCF_YOUNG_GSM1689152.replicated.narrowPeak | sort -k 1,1 -k 2,2n >delete.narrowPeak
bedToBigBed -type=bed6+4 delete.narrowPeak /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes CTCF_replicated_peaks.bigBed
rm delete.narrowPeak


### HIV Integration
awk 'FNR>1 {OFS="\t"; print $2, $3, $3+1, $1, $6, $4, $5, $7}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_integrations.txt | sort -k 1,1 -k 2,2n >delete.bed
bedToBigBed -type=bed6+2 delete.bed /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes HIV_integration.bigBed
rm delete.bed


### HIV Expression; To get a score between 0-1000, had to make negative values positive (sqrt(x^2)), make it an integer (int(x*100))
awk 'FNR>1 {OFS="\t"; print $3, $4, $4+1, $1, sqrt((int($10*100))^2), $5, $2, $6, $7, $8, $9}' /project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.txt | sort -k 1,1 -k 2,2n >delete.bed
bedToBigBed -type=bed6+5 delete.bed /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes HIV_expression.bigBed
rm delete.bed

awk 'FNR>1 {OFS="\t"; print $1, $2, $6}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D.bed | sort -k 1,1 -k 2,2n >delete.bed
bedToBigBed -type=bed3 delete.bed /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes test_loop.bigBed
rm delete.bed

### Hi-ChiP Data; HiChipper
awk '{OFS="\t"}; $8 < 0.05 && $1 == $4 {print $1, $2, $6, "hichipper_rep1", $7, $7, ".", 0, $1, $2, $3, ".", ".", $4, $5, $6, ".", "."}' /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiChipper/YY1_peaks_all/rep1.interactions.all.mango >rep1.bed

awk '{OFS="\t"}; $8 < 0.05 && $1 == $4 {print $1, $2, $6, "hichipper_rep1", $7, $7, ".", 0, $1, $2, $3, ".", ".", $4, $5, $6, ".", "."}' /project/BICF/BICF_Core/shared/Projects/Dorso/Transcription_factor_activity/HiC_pro/HiChipper/YY1_peaks_all/rep2.interactions.all.mango >rep2.bed

cat rep1.bed rep2.bed | sort -k 1,1 -k 2,2n >Hi-ChiP_hichipper.bed

bedToBigBed -as=interact.as -type=bed5+13 Hi-ChiP_hichipper.bed /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes Hi-ChiP_hichipper.bigBed

rm rep1.bed rep2.bed Hi-ChiP_hichipper.bed

### Hi-C files of merged tads and loops; in bed format -> bigBed/bigInteract
awk 'FNR>1 {OFS="\t"}; $1 == $4 {print $1, $2, $6, "HiCloop", "1", "1", ".", 0, $1, $2, $3, ".", ".", $4, $5, $6, ".", "."}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D.bed | sort -k 1,1 -k 2,2n >loop.bed

bedToBigBed -as=interact.as -type=bed5+13 loop.bed /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes Hi-C_loop.bigBed

awk 'FNR>1 {OFS="\t"}; $1 == $4 {print $1, $2, $6, "HiCloop", "1", "1", ".", 0, $1, $2, $3, ".", ".", $4, $5, $6, ".", "."}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.tad.2D.bed | sort -k 1,1 -k 2,2n >tad.bed

bedToBigBed -as=interact.as -type=bed5+13 tad.bed /project/shared/bicf_workflow_ref/human/GRCh38/chrom.sizes Hi-C_tad.bigBed
