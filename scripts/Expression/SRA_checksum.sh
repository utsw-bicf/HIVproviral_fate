#! /bin/sh
### This file runs a test to validate that the files were downloaded correctly

files='/project/BICF/BICF_Core/shared/Projects/Dorso/Expression/rnaseq/raw_data/*.sra'
module load sra_toolkit/2.8.2-1

for file in $files; do
  echo $file >>checksum.txt
  vdb-validate $file >>checksum.txt
done
