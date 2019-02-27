#! /bin/bash
### Process DamID-seq find peaks step
### Used already made pipeline https://github.com/owenjm/find_peaks v 1.0.1 downloaded 2/27/19
#SBATCH --job-name=findpeaks
#SBATCH --partition=super
#SBATCH --output=fp.%j.out
#SBATCH --error=fp.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

perl /work/BICF/s185797/programs/find_peaks-master/find_peaks \
  --fdr=0.05 \
  /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/LaminB11-vs-Dam.gatc.bedgraph
