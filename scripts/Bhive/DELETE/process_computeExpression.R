setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/work/9e/0c85c1024b9ae95676614a3e5eac17')

# Load corplot lib.
source('corplot.R')

# Load integ file.
integs = read.table(file='/project/BICF/BICF_Core/shared/Projects/Dorso/Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_integrations.txt', header=TRUE, comment.char='#')

# Load PCR files.
files = Sys.glob('*.stc')
rna = data.frame(brcd=c(),cnt=c(),rep=c(),tec=c())
dna = data.frame(brcd=c(),cnt=c(),rep=c(),tec=c())
for (f in files) {
  name = strsplit(strsplit(f,'.',fixed=TRUE)[[1]][1],'_')[[1]]
  data = subset(read.table(f,header=FALSE,as.is=TRUE), nchar(V1) < 22 & V2 > 10)
  names(data) = c('brcd','cnt')
  data$rep = as.numeric(name[2])
  data$tec = as.numeric(name[3])
  if (name[1] == 'dna') {
    dna = rbind(dna,data)
  } else if (name[1] == 'rna') { 
    rna = rbind(rna,data)
  }
}

# Merge barcodes by technical replicate.
#??? all = merge(dna,rna,by=c('brcd','rep','tec'))
all = merge(dna,rna,by=c('brcd','rep','tec'), all=T)
all$expr = all$cnt.y/all$cnt.x
all$logexpr = log10(all$expr)

# Iterate over biological replicates
reps = unique(all$rep)
for (r in reps) {
  # Get data from this bio replicate
  repdata = all[all$rep == r,]
  repdna  = dna[dna$rep == r,]
  reprna  = rna[rna$rep == r,]
  
  # Get technical replicate numbers
  tecs = unique(repdata$tec)
  nt = length(tecs)
  
  # Merge all technical replicates and store expr/counts in data frames
  # Expression counts
  dnacnt  = repdna[repdna$tec == tecs[1],][c('brcd','cnt')]
  #rnacnt  = repdna[repdna$tec == tecs[1],][c('brcd','cnt')] ###########!!!!!!!!!! Mistake
  rnacnt  = reprna[reprna$tec == tecs[1],][c('brcd','cnt')]
  # Scatterplot data
  repexpr = repdata[repdata$tec == tecs[1],][c('brcd','logexpr')]
  names(repexpr) = c('brcd',paste('rep',paste(r,tecs[1])))
  
  i = 2
  while (i <= nt) {
    # loop update
    t = tecs[i]
    i = i + 1
    
    # Expression counts
    dnacnt = merge(dnacnt, repdna[repdna$tec == t,][c('brcd','cnt')], by='brcd', all=TRUE)
    rnacnt = merge(rnacnt, reprna[reprna$tec == t,][c('brcd','cnt')], by='brcd', all=TRUE)
    
    # Scatterplot data
    tecdata = repdata[repdata$tec==t,][c('brcd','logexpr')]
    names(tecdata) = c('brcd',paste('rep',paste(r,t,sep='_'),sep='_'))
    repexpr = merge(repexpr,tecdata,by='brcd')
  }
  
  # Expression table by biological replicate
  dnacnt = data.frame(brcd=dnacnt$brcd,count=rowSums(dnacnt[,-1], na.rm=TRUE))
  rnacnt = data.frame(brcd=rnacnt$brcd,count=rowSums(rnacnt[,-1], na.rm=TRUE))
  allcnt = merge(dnacnt,rnacnt,by='brcd')
  names(allcnt) = c('brcd','dna','rna')
 # allcnt$expr = scale(log10(allcnt$rna/allcnt$dna), scale=FALSE)
  allcnt$expr = log10(allcnt$rna/allcnt$dna)
  allcnt$rep  = r
  
  # Merge expression with integration sites
  hivexpr = merge(integs, allcnt, by=c('brcd','rep'))
  
  # Write table
  write.table(hivexpr, file=("testing_hiv_expression.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, append=TRUE)

  # Create figures dir
#  dir.create('figures', showWarnings = FALSE)
  
  # Cor plot
 # pdf(paste(paste('figures/expr_rep',r,sep='_'),'pdf',sep='.'))
 # corplot(repexpr[,-1])
 # dev.off()
}
