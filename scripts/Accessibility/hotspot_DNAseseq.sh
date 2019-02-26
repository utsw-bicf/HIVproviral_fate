#! /bin/bash
### This script runs hotspot2 to identify if the DNAse-seq data is quality

module load bedops/2.4.14
module load samtools/1.6
module load UCSC_userApps/v317
module load gcc/6.1.0

### First run scripts/extractCenterSites.sh, note that in the chrom.sizes file the second column needs to be a 0
### Blacklist was downloaded from http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/
/work/BICF/s185797/programs/hotspot2-master/scripts/extractCenterSites.sh -c chrom_sizes.bed -o ecs.starch -M /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis/workflow/output_DNAse/hotspot2/hg38_whitelist.bed

### Run hotspot2
/work/BICF/s185797/programs/hotspot2-master/scripts/hotspot2.sh -c chrom_sizes.bed -C ecs.starch /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis/workflow/output_DNAse/filterReads/ENCFF001DPF.filt.nodup.bam ENCFF001DPF

/work/BICF/s185797/programs/hotspot2-master/scripts/hotspot2.sh -c chrom_sizes.bed -C ecs.starch /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis/workflow/output_DNAse/filterReads/ENCFF001DPG.filt.nodup.bam ENCFF001DPG

### the *SPOT.txt contains the spot score, which should be >0.4

