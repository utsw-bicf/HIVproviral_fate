#! /bin/sh
### Download RNAseq data for d'Orso lab
### Chipseq

#SBATCH --job-name=dd
#SBATCH --partition=super
#SBATCH --output=dd.%j.out
#SBATCH --error=dd.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=END

module load sra_toolkit/2.8.2-1

#SRAs="SRR061744 SRR074203 SRR1057276 SRR1522115 SRR1522118 SRR1603646 SRR1603648 SRR1603649 SRR1603650 SRR1603652 SRR2043612 SRR2043614 SRR3722570 SRR3722573 SRR577482 SRR577483 SRR6010201 SRR6010202"
#SRAs="SRR3722566"
SRAs="SRR1791460 SRR1791461 SRR1791462 SRR1791463 SRR1791464 SRR1791465"

for sra in $SRAs; do
  sixsra=`echo ${sra} | cut -c1-6`
  wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${sixsra}/${sra}/${sra}.sra
done

