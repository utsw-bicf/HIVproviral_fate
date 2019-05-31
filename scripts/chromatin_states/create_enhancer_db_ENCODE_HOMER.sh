#! /bin/bash
### Use encode data to identify superenhancers

### CD4
#CL:0000897\|CL:0000624\|CL:0000792\|CL:0000905\|CL:0000895

### T-cell
#CL:0000084
#CL:0000899

########## Parse the database from Chris
grep "CL:0000897\|CL:0000624\|CL:0000792\|CL:0000905\|CL:0000895" /project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/get_encode_run_only_ChIP/GRCh38/metadata_GRCh38_20190523.tsv | grep "GRCh38" | grep -i "chip" | grep -i "peaks" | awk '{OFS="\t"} $2 != "bam" {print $0}' | grep -v "extremely low read depth" | awk -F "\t" '{OFS="\t"} $3 = "replicated peaks" {print $0}' | awk -F "\t" '{OFS="\t"} $31 ~ /2/ {print $0}' >/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/enhancer_database/ChrisBennet_ENCODE/get_encode_run_only_ChIP/GRCh38/test.tsv

########## Replicate peaks of the replicate peak files

