#! /bin/sh
### This checks that the ENCODE data was downloaded properly

echo "ENCFF001DPG" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/raw_data/checksum.txt
echo "MD5sum	e026f2105135789d8f2262b68a533230" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/raw_data/checksum.txt
echo "Content MD5sum	c36f4b30164be147c54c1676e5f76fda" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/raw_data/checksum.txt

echo "ENCFF001DPF" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/raw_data/checksum.txt
echo "MD5sum	a3ccba847ee9b8e7467bf0a60151af38" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/raw_data/checksum.txt
echo "Content MD5sum	62b975ef1ccae2b8236aaac6f0ef4f1b" >>/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/raw_data/checksum.txt

md5sum /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/raw_data/ENCFF001DPG.fastq.gz >>/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/raw_data/checksum.txt
md5sum /project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/raw_data/ENCFF001DPF.fastq.gz >>/project/BICF/BICF_Core/shared/Projects/Dorso/Accessibility/raw_data/checksum.txt
