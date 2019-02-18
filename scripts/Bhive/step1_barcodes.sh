#! /bin/sh

###Programs used
# seeq-1.1 (commit a6b091b28751b8dc7ec9ed2563e82f2aeb5207c7, downloaded 2/15/19 onto local mac)
# starcode-v1.3 17-07-2018 (commit 11c0b0f12167c20e84d1a97bd2185d3478f4dde5, downloaded 2/17/19 on local mac)
### Detect the whole line and barcodes for integration sites
#gunzip -c /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/SRR3614675_1.fastq.gz | /Users/holly/Desktop/seeq/seeq/seeq -d 3 TATAGTGAGTCGTA >/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_1_wholeline.txt
#gunzip -c /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/SRR3614675_1.fastq.gz | /Users/holly/Desktop/seeq/seeq/seeq -d 3 -r TATAGTGAGTCGTA >/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_1_barcode.txt

#gunzip -c /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/SRR3614676_1.fastq.gz | /Users/holly/Desktop/seeq/seeq/seeq -d 3 TATAGTGAGTCGTA >/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_1_wholeline.txt
#gunzip -c /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/SRR3614676_1.fastq.gz | /Users/holly/Desktop/seeq/seeq/seeq -d 3 -r TATAGTGAGTCGTA >/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_1_barcode.txt

### Detect errors in barcodes
#/Users/holly/Desktop/starcode/starcode/starcode -d 1 --print-clusters -i /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_1_barcode.txt -o /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_1_starcode.txt
#/Users/holly/Desktop/starcode/starcode/starcode -d 1 --print-clusters -i /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_1_barcode.txt -o /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_1_starcode.txt

### match the pulled barcodes with the reverse (_2) fastq sequence
# Remove duplicates from whole file
#sort /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_1_wholeline.txt | uniq >/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_1_wholeline_unique.txt
#sort /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_1_wholeline.txt | uniq >/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_1_wholeline_unique.txt

# Pull the sequence header from the matches
while read LINE; do
    gunzip -c /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/SRR3614675_1.fastq.gz | grep -B 1 -A 2 "$LINE" | grep -v -e "--" >>/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_1_keep.fastq
done </Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_1_wholeline_unique.txt
gzip /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_1_keep.fastq >/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_1_keep.fastq.gz
rm /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_1_keep.fastq

while read LINE2; do
    gunzip -c /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/SRR3614676_1.fastq.gz | grep -B 1 -A 2 "$LINE2" | grep -v -e "--" >>/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_1_keep.fastq
done </Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_1_wholeline_unique.txt
gzip /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_1_keep.fastq >/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_1_keep.fastq.gz
rm /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_1_keep.fastq

# Pull the header line from fastq1
gunzip -c /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_1_keep.fastq | grep "^@N" >/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_1_header.txt
gunzip -c /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_1_keep.fastq | grep "^@N" >/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_1_header.txt

# Pull fastq2 sequence from header list
while read HDLINE; do
    gunzip -c /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/SRR3614675_2.fastq.gz | grep -A 3 "$HDLINE" | grep -v -e "--" >>/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_2_keep.fastq.gz
done </Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_1_header.txt
gzip /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_2_keep.fastq >/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_2_keep.fastq.gz
rm /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614675_2_keep.fastq

while read HDLINE2; do
    gunzip -c /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/raw_data/SRR3614676_2.fastq.gz | grep -A 3 "$HDLINE2" | grep -v -e "--" >>/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_2_keep.fastq
done </Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_1_header.txt
gzip /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_2_keep.fastq >/Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_2_keep.fastq.gz
rm /Volumes/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/barcodes/SRR3614676_2_keep.fastq

