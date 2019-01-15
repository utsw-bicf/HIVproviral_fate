#! /bin/sh
### Download RNAseq data for d'Orso lab
### Chipseq

#SBATCH --job-name=dd
#SBATCH --partition=super
#SBATCH --output=dd.%j.out
#SBATCH --error=dd.%j.err
#SBATCH --mail-user=holly.ruess@utsouthwestern.edu
#SBATCH --mail-type=NONE

module load sra_toolkit/2.8.2-1

### GSM569085: SRR061743
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR061/SRR061743/SRR061743.sra

### GSM945268: SRR577484
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR577/SRR577484/SRR577484.sra

### GSM1519642: SRR1603654
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR160/SRR1603654/SRR1603654.sra

### GSM569089: SRR061747
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR061/SRR061747/SRR061747.sra

### GSM1062740: SRR647929
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR647/SRR647929/SRR647929.sra

### GSM1850204: SRR2157609
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR215/SRR2157609/SRR2157609.sra

### GSM2218766: SRR3722577
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR372/SRR3722577/SRR3722577.sra

### GSM1603229: SRR1791514
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791514/SRR1791514.sra

### GSM1603229: SRR1791515
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791515/SRR1791515.sra

### GSM1603229: SRR1791516
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791516/SRR1791516.sra

### GSM1603229: SRR1791517
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791517/SRR1791517.sra

### GSM1603229: SRR1791518
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791518/SRR1791518.sra

### GSM1603229: SRR1791519
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791519/SRR1791519.sra

### GSM1603225: SRR1791493
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791493/SRR1791493.sra

### GSM1603225: SRR1791494
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791494/SRR1791494.sra

### GSM1603225: SRR1791495
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791495/SRR1791495.sra

### GSM1603225: SRR1791496
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791496/SRR1791496.sra

### GSM1603225: SRR1791497
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791497/SRR1791497.sra

### GSM1603225: SRR1791498
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791498/SRR1791498.sra

### GSM1603213: SRR1791425
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791425/SRR1791425.sra

### GSM1603213: SRR1791426
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791426/SRR1791426.sra

### GSM1603213: SRR1791427
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791427/SRR1791427.sra

### GSM1603213: SRR1791428
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791428/SRR1791428.sra

### GSM1603213: SRR1791429
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791429/SRR1791429.sra

### GSM1603213: SRR1791430
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791430/SRR1791430.sra

### GSM1603209: SRR1791404
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791404/SRR1791404.sra

### GSM1603209: SRR1791405
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791405/SRR1791405.sra

### GSM1603209: SRR1791406
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791406/SRR1791406.sra

### GSM1603209: SRR1791407
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791407/SRR1791407.sra

### GSM1603209: SRR1791408
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791408/SRR1791408.sra

### GSM1603209: SRR1791409
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791409/SRR1791409.sra

### GSM1603215: SRR1791437
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791437/SRR1791437.sra

### GSM1603215: SRR1791438
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791438/SRR1791438.sra

### GSM1603215: SRR1791439
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791439/SRR1791439.sra

### GSM1603215: SRR1791440
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791440/SRR1791440.sra

### GSM1603215: SRR1791441
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791441/SRR1791441.sra

### GSM1603215: SRR1791442
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791442/SRR1791442.sra

### GSM1603227: SRR1791505
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791505/SRR1791505.sra

### GSM1603227: SRR1791506
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791506/SRR1791506.sra

### GSM1603227: SRR1791507
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791507/SRR1791507.sra

### GSM1603227: SRR1791508
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791508/SRR1791508.sra

### GSM1603227: SRR1791509
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791509/SRR1791509.sra

### GSM1603227: SRR1791510
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791510/SRR1791510.sra

### GSM1603211: SRR1791413
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791413/SRR1791413.sra

### GSM1603211: SRR1791414
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791414/SRR1791414.sra

### GSM1603211: SRR1791415
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791415/SRR1791415.sra

### GSM1603211: SRR1791416
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791416/SRR1791416.sra

### GSM1603211: SRR1791417
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791417/SRR1791417.sra

### GSM1603211: SRR1791418
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791418/SRR1791418.sra

### GSM628535: SRR074201
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR074/SRR074201/SRR074201.sra

### GSM628532: SRR074195
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR074/SRR074195/SRR074195.sra

### GSM628532: SRR074196
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR074/SRR074196/SRR074196.sra

### GSM1603217: SRR1791449
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791449/SRR1791449.sra

### GSM1603217: SRR1791450
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791450/SRR1791450.sra

### GSM1603217: SRR1791451
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791451/SRR1791451.sra

### GSM1603217: SRR1791452
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791452/SRR1791452.sra

### GSM1603217: SRR1791453
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791453/SRR1791453.sra

### GSM1603217: SRR1791454
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791454/SRR1791454.sra

### GSM1603223: SRR1791484
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791484/SRR1791484.sra

### GSM1603223: SRR1791485
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791485/SRR1791485.sra

### GSM1603223: SRR1791486
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791486/SRR1791486.sra

### GSM1603223: SRR1791487
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791487/SRR1791487.sra

### GSM1603223: SRR1791488
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791488/SRR1791488.sra

### GSM1603223: SRR1791489
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791489/SRR1791489.sra

### GSM1603221: SRR1791472
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791472/SRR1791472.sra

### GSM1603221: SRR1791473
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791473/SRR1791473.sra

### GSM1603221: SRR1791474
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791474/SRR1791474.sra

### GSM1603221: SRR1791475
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791475/SRR1791475.sra

### GSM1603221: SRR1791476
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791476/SRR1791476.sra

### GSM1603221: SRR1791477
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1791477/SRR1791477.sra

### GSM2218755: SRR3722566
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR3722566/SRR3722566.sra

### GSM1442004: SRR1522115
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1522115/SRR1522115.sra

### GSM2773998: SRR6010201
### single-end
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR6010201/SRR6010201.sra


###
### fastq-dump --gzip --skip-technical  --readids --dumpbase --split-files --clip SRR6010201

###Recommended by dorso lab
### fastq-dump --gzip --outdir path/to/dir --dumpbase --origfmt input.fi

