#! /bin/sh
### Run Hi-C data through HOMER

#SBATCH --job-name=stdHOMER
#SBATCH -p 256GB 
#SBATCH --mem 253952
#SBATCH --output=stdHOMER.%j.out
#SBATCH --error=stdHOMER.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### Follow: http://homer.ucsd.edu/homer/interactions2/quickndirty.html
#module load homer/4.10.4  #### Was installed locally on /project/BICF/BICF_Core/s185797/homer; and in bashrc
module load juicebox/1.5.6
module load R

### Concatenate fastqs by library
echo "Start concatenating libraries"
## lib1
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8244644_1.fastq.gz \
/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8244645_1.fastq.gz \
>/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489136_1.fastq.gz

cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8244644_2.fastq.gz \
/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8244645_2.fastq.gz \
>/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489136_2.fastq.gz

## lib2
cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8244646_1.fastq.gz \
/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388751_1.fastq.gz \
#project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388752_1.fastq.gz \
/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388753_1.fastq.gz \
/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388754_1.fastq.gz \
/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388755_1.fastq.gz \
/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388756_1.fastq.gz \
/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388757_1.fastq.gz \
>/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489137_1.fastq.gz

cat /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8244646_2.fastq.gz \
/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388751_2.fastq.gz \
/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388752_2.fastq.gz \
/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388753_2.fastq.gz \
/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388754_2.fastq.gz \
/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388755_2.fastq.gz \
/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388756_2.fastq.gz \
/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/SRR8388757_2.fastq.gz \
>/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489137_2.fastq.gz
echo "End concatenating libraries"

### Trim reads
echo "Start read trimming"
fastqs='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489136_1.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489136_2.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489137_1.fastq.gz /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489137_2.fastq.gz'

for fastq in ${fastqs}; do
  homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 ${fastq}
done
echo "End read trimming"

### Make reference genome for bowtie using only canonical chromosomes (1-22,X,Y,M)
echo "Start making canonical reference"
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/canonical_genome

module load samtools
samtools faidx /project/shared/bicf_workflow_ref/human/GRCh38/genome.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/canonical_genome/canonical_genome.fa
module rm samtools

module load bowtie2/2.2.8-intel
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/canonical_genome/bowtie2_index/
bowtie2-build -f /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/canonical_genome/canonical_genome.fa /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/canonical_genome/bowtie2_index/genome
echo "End making canonical reference"

### Alignment with bowtie2
echo "Start alignment"
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/alignments
tfastqs='/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489136_1.fastq.gz.trimmed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489136_2.fastq.gz.trimmed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489137_1.fastq.gz.trimmed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/raw_data/GSM3489137_2.fastq.gz.trimmed'

for tfastq in ${tfastqs}; do
  outname=$(basename ${tfastq} .fastq.gz.trimmed)
  bowtie2 -p 20 -x /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/canonical_genome/bowtie2_index/genome -U ${tfastq} >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/alignments/${outname}.sam
done
module rm bowtie2/2.2.8-intel
echo "End alignment"

### HI-C tag directory
echo "Start Hi-C tagging"
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2
makeTagDirectory std_HicHomerTagDir/lib1 alignments/GSM3489136_1.sam,alignments/GSM3489136_2.sam -genome hg38 -checkGC -restrictionSite GATC 
makeTagDirectory std_HicHomerTagDir/lib2 alignments/GSM3489137_1.sam,alignments/GSM3489137_2.sam -genome hg38 -checkGC -restrictionSite GATC 
echo "End Hi-C tagging"

### Create JuiceBox *.hic
### Uses JuiceBox
### If the file is big use the export function for perl
echo "Create JuiceBox file"
tagDir2hicFile.pl \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1 \
  -juicer auto \
  -genome hg38 \
  -p 20 \
  -juicerExe "java -jar /cm/shared/apps/juicebox/1.5.6/juicer_tools_linux_0.8.jar" 
export TMPDIR=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/tmp
tagDir2hicFile.pl \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2 \
  -juicer auto \
  -genome hg38 \
  -p 20 \
  -juicerExe "java -jar /cm/shared/apps/juicebox/1.5.6/juicer_tools_linux_0.8.jar"
echo "End create JuiceBox file"

##### Visualize a Hi-C contact map;   -bgonly \
echo "Create Background Models"
analyzeHiC \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1 \
  -cpu 8
analyzeHiC \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2 \
  -cpu 8
echo "End Create Background Models"


### Compartment analysis
### May need to do this again after runnning chromHMM to find open and closed regions
echo "Compartment analsysis - PCA"
export TMPDIR=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/tmp
runHiCpca.pl auto \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1 \
  -cpu 10
export TMPDIR=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/tmp
runHiCpca.pl auto \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2 \
  -cpu 10
###combine outputs
echo "combine PCA outputs"
annotatePeaks.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2/lib2.50x100kb.PC1.txt \
  hg38 \
  -noblanks \
  -bedGraph /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2/lib2.50x100kb.PC1.bedGraph /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1/lib1.50x100kb.PC1.bedGraph \
  > /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/combine_PC1_output.txt
echo "Compartment analsysis - PCA"

### Chromatin Compaction
echo "Chromatin Compaction"
analyzeHiC \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1 \
  -res 5000 \
  -window 15000 \
  -nomatrix \
  -compactionStats auto \
  -cpu 10
analyzeHiC \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2 \
  -res 5000 \
  -window 15000 \
  -nomatrix \
  -compactionStats auto \
  -cpu 10
echo "End Chromatin Compaction"


### Finding TADs and Loops
### To get bad regions (done on 4/19/19):
wget -O gap.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz
wget -O dups.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz
zcat gap.txt.gz dups.txt.gz | cut -f2-4 > badRegions.bed
echo "Start finding TADs and Loops"
export TMPDIR=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/tmp
findTADsAndLoops.pl \
  find /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1 \
  -cpu 10 \
  -res 3000 \
  -window 15000 \
  -genome hg38 \
  -p /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/badRegions.bed

export TMPDIR=/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/tmp
findTADsAndLoops.pl \
  find /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2 \
  -cpu 10 \
  -res 3000 \
  -window 15000 \
  -genome hg38 \
  -p /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/badRegions.bed
echo "End finding TADs and Loops"

####################################################################################################

### Make interaction matrix
 module load juicebox/1.5.6
 java -jar /cm/shared/apps/juicebox/1.5.6/Juicebox.jar
#load *.hic files

### Merge tads and loops
merge2Dbed.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1/lib1.tad.2D.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2/lib2.tad.2D.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.tad.2D.bed

merge2Dbed.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1/lib1.loop.2D.bed /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2/lib2.loop.2D.bed >/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D.bed

### Find Significant Interactions
analyzeHiC /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1/ -res 1000000 -interactions lib1_significantInteractions.txt -nomatrix -cpu 25
analyzeHiC /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2/ -res 1000000 -interactions lib2_significantInteractions.txt -nomatrix -cpu 25

### Find High resolution connections
findHiCInteractionsByChr.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib1/ -res 2000 -superRes 10000 -cpu 25 > lib1_chrInteractions.txt
findHiCInteractionsByChr.pl /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/lib2/ -res 2000 -superRes 10000 -cpu 25 > lib2_chrInteractions.txt

### Annotate interactions
annotateInteractions.pl lib1_chrInteractions.txt hg38 Annotated_interations_lib1/ -pvalue 0.05 
#annotateInteractions.pl lib2_chrInteractions.txt hg38 Annotated_interations_lib2/ -pvalue 0.05 

### Find interactions with peaks for cytoscape
annotateInteractions.pl lib1_chrInteractions.txt hg38 lib1_annotated_w_peaks/ -p /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/consensusPeaks/H3K4me3.replicated.narrowPeak /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/consensusPeaks/H3K27ac.replicated.narrowPeak /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/GSM1603225_GFP_H3K4me1.replicated.narrowPeak /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/GSM1603209_GFP_H3K36me3.replicated.narrowPeak /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/GSM1603215_GFP_H3K79me3.replicated.narrowPeak /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/GSM1603227_GFP_H3K9me3.replicated.narrowPeak

annotateInteractions.pl lib2_chrInteractions.txt hg38 lib2_annotated_w_peaks/ -p /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/consensusPeaks/H3K4me3.replicated.narrowPeak /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_consensus/consensusPeaks/H3K27ac.replicated.narrowPeak /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/GSM1603225_GFP_H3K4me1.replicated.narrowPeak /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/GSM1603209_GFP_H3K36me3.replicated.narrowPeak /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/GSM1603215_GFP_H3K79me3.replicated.narrowPeak /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/consensusPeaks/GSM1603227_GFP_H3K9me3.replicated.narrowPeak
