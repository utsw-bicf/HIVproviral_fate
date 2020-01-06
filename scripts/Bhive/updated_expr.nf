def print_help = {
    log.info ''
    log.info '     BHIVE expression pipeline'
    log.info '-----------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    expression.nf --integs HIV_INTEGS [optional arguments]'
    log.info ''
    log.info 'Computes HIV expression from PCR and RT-PCR reads.'
    log.info 'Additional params, such as datasets for each replicate,'
    log.info 'flanking sequences, etc. should be defined in "params.txt".'
    log.info ''
    log.info 'Required arguments:'
    log.info ''
    log.info '  --integs      Path to file containing HIV integrations.'
    log.info ''
    log.info 'Optional arguments:'
    log.info '  --options     Path to alternative mapping params file (default: params.txt)'
    log.info '  --out         Output path (default: .)'
    log.info ''
}

lib_dir = "lib"

// Set defaults
params.help        = false
params.out         = '.'
params.options     = 'expr.cfg'
params.integs      = null

/**********************************************************************
**************************** PARSE  PARAMS ****************************
**********************************************************************/

if (params.help) {
   print_help()
   exit 0
}

if (params.integs == null) {
   log.info "error: HIV integration file not specified!"
   exit 1
}

ifile = file("${params.integs}")
if (!ifile.isFile()) {
   log.info "error: HIV integration file not found in '${ifile}'."
   exit 1
}
integ_file = Channel.from(ifile)

/**********************************************************************
***************************** PARSE OPTIONS ***************************
**********************************************************************/
/*
** Mapping options format:
** 
** 
** datasets:
** biological replicate,technical replicate,filename/URL/{SRR,ERR,DRR} reference
** 
*/

def options = ['bfs': null]
def names = ['','dna','rna']

mfile = file("${params.options}")
if (!mfile.isFile()) {
   log.info "error: options file not found (${mfile})!"
   exit 1
}

// Varaibles
def status = 0

// Channels
file_getref = Channel.create()
datasets    = Channel.create()

// Parse file by lines.
mfile.eachLine { line ->
   if (line =~ /^#/ || line =~ /^\s*$/) {
      return null
   } else if (line =~ /^dna:/) {
      status = 1
      return null
   } else if (line =~ /^rna:/) {
      status = 2
      return null
   }
   switch (status) {
   case 0:
     t = line.split(/\s*=\s*/).toList()
     if (!options.containsKey(t[0])) {
        log.info "unknown option: '$t[0]' in $params.options"
        exit 1
     }
     options[t[0]] = t[1]
     break
   case 1..2:
     t = line.split(/\s*,\s*/).toList()
     if (t.size() == 3) {
        if (t[-1] =~ /^SRR|^ERR|^DRR/ || t[-1] =~ /^http|^ftp/) {
           file_getref << [ t[-1], [names[status],t[0],t[1]] ]
        } else {
           reads = file("${t[-1]}")
           if (reads.isFile()) {
              datasets << [ reads, [names[status],t[0],t[1]] ]
           } else {
              log.info "error: file not found, in '$line'. (URLs must start with http/ftp)"
              exit 1
           }
        }
     } else {
        log.info "error in entry '$line': rna/dna entries must specify a file path, an URL or a GEO reference starting with SRR/ERR/DRR."
        exit 1
     }
     break
   }
}

// Parse options
if (!options['bfs']) {
   log.info "error: 'bfs' option must be defined in $params.options before 'dna:' or 'rna:'."
   exit 1
}

// Group the same references/URLS to avoid redundant downloads.
file_getref.close()

process getDataset {
   // Process options
   tag "${data[0]}_${data[1]}"
   publishDir path:"${params.out}/datasets/", mode:'symlink'
   // Cluster options
   cpus 1
   memory '2GB'

   input:
   val data from file_getref
   output:
   set file('*.fastq.gz'), sample into datasets
   script:
   sample = data[1]
   ref = data[0]
   if (ref =~ /^SRR|^ERR|^DRR/) {
      """
      fastq-dump --split-files --gzip -A ${ref}
      rm -f ~/ncbi/public/sra/${ref}.sra
      """
   } else {
      """
      wget ${ref}
      """
   }
}

/**********************************************************************
************************* EXPRESSION PIPELINE *************************
**********************************************************************/


process extractPCRBarcodes {
   // Process options
   tag "${sample[0]}_${sample[1]}_${sample[2]}"
   publishDir path:"${params.out}/barcodes_raw/", mode:'symlink'
   
   // Custer options
   cpus 1
   memory '4GB'

   input:
   set file(reads), sample from datasets

   output:
   set file("*.bcd"), sample into pcr_barcodes
   
   script:
   outf = "${sample[0]}_${sample[1]}_${sample[2]}.bcd"
"""
#!/usr/bin/env python
import seeq
import gzip
MIN_BRCD_LEN = 15

# Helper functions.
def getbarcode(matcher, txt):
   return matcher.matchPrefix(txt, False) or ''

if __name__ == '__main__':
   fname = '${reads}'
   fout = '${outf}'
   outp = open(fout, 'w+')
   T7seq = '${options['bfs']}'
   T7 = seeq.compile(T7seq, int(max(1,round(0.2*len(T7seq)))))
   with gzip.open(fname) as f:
      # Read fastq file.
      for lineno,line in enumerate(f):
         if lineno % 4 != 1: continue
         barcode = getbarcode(T7, line)
         if len(barcode) < MIN_BRCD_LEN: continue
         outp.write(barcode + '\\n')

   outp.close()

"""
}

process clusterPCRBarcodes {
   // Process options
   tag "$sample[1]"
   publishDir path:"${params.out}/barcodes_cluster/", mode:'symlink'
   
   // Cluster options
   cpus 12
   memory '64GB'

   input:
   set file(barcodes), sample from pcr_barcodes
   
   output:
   file "*.stc" into pcr_clusters
   
   script:
   """
   starcode -d ${options['dist']} -qt ${task.cpus} ${barcodes} > ${sample[0]}_${sample[1]}_${sample[2]}.stc
   """
}

process computeExpression {
   // Process options
   tag ""
   publishDir path:"${params.out}/", mode:'move'
   // Cluster options
   cpus 1
   memory '4GB'

   input:
   file barcodes from pcr_clusters.flatten().toList()
   file integ_path from integ_file
   file libs from Channel.fromPath("${lib_dir}/*.R").flatten().toList()
   //   file hgene_path from annotation_file_exp
   //   file rnaseq_path from rnaseq_data

   output:
   file 'hiv_expression.txt' into hiv_expression
   file 'figures/*.pdf' into expression_figures

   script:
   """
   #!/usr/bin/env Rscript

   # Load corplot lib.
   source('corplot.R')

   # Load integ file.
   integs = read.table(file='${integ_path}', header=TRUE, comment.char='#')

   # Load PCR files.
   files = Sys.glob('*.stc')
   rna = data.frame(brcd=c(),cnt=c(),rep=c(),tec=c())
   dna = data.frame(brcd=c(),cnt=c(),rep=c(),tec=c())
   for (f in files) {
     name = strsplit(strsplit(f,'.',fixed=TRUE)[[1]][1],'_')[[1]]
     data = subset(read.table(f,header=FALSE,as.is=TRUE), nchar(V1) < 22 & V2 > 10)
     names(data) = c('brcd','cnt')
     data\$rep = as.numeric(name[2])
     data\$tec = as.numeric(name[3])
     if (name[1] == 'dna') {
       dna = rbind(dna,data)
     } else if (name[1] == 'rna') { 
       rna = rbind(rna,data)
     }
   }

   # Merge barcodes by technical replicate.
   all = merge(dna,rna,by=c('brcd','rep','tec'))
   all\$expr = all\$cnt.y/all\$cnt.x
   all\$logexpr = log10(all\$expr)

   # Iterate over biological replicates
   reps = unique(all\$rep)
   for (r in reps) {
     # Get data from this bio replicate
     repdata = all[all\$rep == r,]
     repdna  = dna[dna\$rep == r,]
     reprna  = rna[rna\$rep == r,]

     # Get technical replicate numbers
     tecs = unique(repdata\$tec)
     nt = length(tecs)

     # Merge all technical replicates and store expr/counts in data frames
     # Expression counts
     dnacnt  = repdna[repdna\$tec == tecs[1],][c('brcd','cnt')]
     rnacnt  = reprna[reprna\$tec == tecs[1],][c('brcd','cnt')]
     # Scatterplot data
     repexpr = repdata[repdata\$tec == tecs[1],][c('brcd','logexpr')]
     names(repexpr) = c('brcd',paste('rep',paste(r,tecs[1])))

     i = 2
     while (i <= nt) {
       # loop update
       t = tecs[i]
       i = i + 1

       # Expression counts
       dnacnt = merge(dnacnt, repdna[repdna\$tec == t,][c('brcd','cnt')], by='brcd', all=TRUE)
       rnacnt = merge(rnacnt, reprna[reprna\$tec == t,][c('brcd','cnt')], by='brcd', all=TRUE)

       # Scatterplot data
       tecdata = repdata[repdata\$tec==t,][c('brcd','logexpr')]
       names(tecdata) = c('brcd',paste('rep',paste(r,t,sep='_'),sep='_'))
       repexpr = merge(repexpr,tecdata,by='brcd')
     }

     # Expression table by biological replicate
     dnacnt = data.frame(brcd=dnacnt\$brcd,count=rowSums(dnacnt[,-1], na.rm=TRUE))
     rnacnt = data.frame(brcd=rnacnt\$brcd,count=rowSums(rnacnt[,-1], na.rm=TRUE))
     allcnt = merge(dnacnt,rnacnt,by='brcd')
     names(allcnt) = c('brcd','dna','rna')
     allcnt\$expr = scale(log10(allcnt\$rna/allcnt\$dna), scale=FALSE)
     allcnt\$rep  = r

     # Merge expression with integration sites
     hivexpr = merge(integs, allcnt, by=c('brcd','rep'))

     # Write table
     write.table(hivexpr, file='hiv_expression.txt', sep='\\t', quote=FALSE, row.names=FALSE, col.names=TRUE, append=TRUE)

     # Create figures dir
     dir.create('figures', showWarnings = FALSE)

     # Cor plot
     pdf(paste(paste('figures/expr_rep',r,sep='_'),'pdf',sep='.'))
     corplot(repexpr[,-1])
     dev.off()
   }
   """
}
