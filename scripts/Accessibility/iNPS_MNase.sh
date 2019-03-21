#!/bin/bash

#SBATCH --job-name=iNPS
#SBATCH --partition=super
#SBATCH --output=iNPS.%j.out
#SBATCH --error=iNPS.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

### This script is an MNase-Seq caller using iNPS (http://www.picb.ac.cn/hanlab/iNPS.html)
### iNPS was downloaded and installed locally in /work/s185797/programs/iNPS_V1.2.2.py

########## Turn bam into 3 column bed; where col 2 and 3 are end to end of the paired reads; keep only chr 1-22, X, Y, M
###### Turn tagAlign file into 3 column bed
zcat /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/convertReads/Jkt106_MNase_Seq_20170218_1_val_1Jkt106_MNase_Seq_20170218_2_val_2.filt.nodup.tn5.tagAlign.gz >/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/iNPS/Jkt106_MNase_Seq_20170218_1_val_1Jkt106_MNase_Seq_20170218_2_val_2.filt.nodup.tn5.tagAlign.txt
paste -d " " - - </project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/iNPS/Jkt106_MNase_Seq_20170218_1_val_1Jkt106_MNase_Seq_20170218_2_val_2.filt.nodup.tn5.tagAlign.txt >/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/iNPS/Jkt106_MNase_seq.txt

awk 'BEGIN {OFS="\t"}; {if ($1==$7 && $2>$8) {print $1,$8,$2;} else if ($1==$7 && $2<$8) print $1,$2,$8;}' /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/iNPS/Jkt106_MNase_seq.txt | grep -v "_" | grep -v "HLA" | sort -k 1,2 -V -s >/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/iNPS/Jkt106_MNase_seq.bed

##### Run iNPS
module load python/3.4.x-anaconda
mkdir /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/iNPS_output
python /work/BICF/s185797/programs/iNPS_V1.2.2.py -i /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/iNPS/Jkt106_MNase_seq.bed --s_p p -o /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/iNPS_output

rm /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/iNPS/Jkt106_MNase_Seq_20170218_1_val_1Jkt106_MNase_Seq_20170218_2_val_2.filt.nodup.tn5.tagAlign.txt /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/atacseq_analysis_MNase/workflow/output_MNase/iNPS/Jkt106_MNase_seq.txt

