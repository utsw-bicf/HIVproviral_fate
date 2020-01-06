#! /bin/sh
### This script combines the fingerprint plots into 1 image

montage \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/*_fingerprint.png \
  -tile 2x13 \
  -mode concatenate \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Epigenetics/chipseq_analysis/workflow/output_011519/experimentQC/fingerprint_plots.pdf
