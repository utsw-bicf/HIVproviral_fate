public boolean isInteger( String input ) {
    try {
        Integer.parseInt( input );
        return true;
    }
    catch( NumberFormatException e ) {
        return false;
    }
}

def print_help = {
    log.info ''
    log.info '     BHIVE mapping pipeline'
    log.info '--------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    mapping.nf --index BWA_INDEX [optional arguments]'
    log.info ''
    log.info 'Computes HIV integration sites from iPCR NGS reads.'
    log.info 'Additional params, such as datasets for each replicate, matching distances,'
    log.info 'mapping qualities, etc. should be defined in "params.txt".'
    log.info ''
    log.info 'Required arguments:'
    log.info ''
    log.info '  --index       Path to BWA index file (required)'
    log.info ''
    log.info 'Optional arguments:'
    log.info '  --options     Path to alternative mapping params file (default: params.txt)'
    log.info '  --out         Output path (default: .)'
    log.info ''
}

script_dir = "src"

extract_barcode_src = Channel.fromPath("${script_dir}/mapping_extract_barcodes.py")
barcode_seqnames_src = Channel.fromPath("${script_dir}/mapping_barcode_seqnames.py")
best_assignment_src = Channel.fromPath("${script_dir}/mapping_best_assignment.py")


// Set defaults
params.help        = false
params.out         = '.'
params.options     = 'map.cfg'
params.index       = null

if (params.help) {
   print_help()
   exit 0
}

/**********************************************************************
**************************** PARSE  PARAMS ****************************
**********************************************************************/


bwa_index   = Channel.create()
index_files = Channel.create()
// 0. Find BWA index
if (params.index) {
   // Check bwa index files.
   fasta_ref = file("${params.index}.bwt")
   index_ref = file("${params.index}")
   if (fasta_ref.isFile()) {
      if (!index_ref.isFile()) {
         index_ref.toFile().createNewFile()
      }
      bwa_index << file("${params.index}")
      index_files << file("${index_ref}.{bwt,amb,ann,pac,sa}")
   } else if (index_ref.getExtension() in ['bwt','amb','ann','pac','sa'] && index_ref.isFile()) {
      index_str = (index_ref.getParent().equals(null) ? './' : index_ref.getParent()) + index_ref.getBaseName()
      index_ref = file(index_str)
      if (!index_ref.isFile()) {
         index_ref.toFile().createNewFile()
      }
      bwa_index << index_ref
      index_files << file("${index_ref}.{bwt,amb,ann,pac,sa}")
   } else {
      log.info "error: BWA index not found in '${params.index}'."
      exit 1
   }
} else {
   log.info "error: '--index' option not specified."
   exit 1
}
bwa_index.close()
index_files.close()


/**********************************************************************
***************************** PARSE OPTIONS ***************************
**********************************************************************/
/*
** Mapping options format:
** 
**  BFS (Barcode flanking sequence)
**  LTR (HIV LTR sequence, must be present in reverse read)
**  RFS (Restriction enzyme flanking sequence, in HIV construct)
**  DIST (Barcode clustering distance)
**  MAPQ (Minimum mapping quality, 20)
**  INTQ (Minimum assignment score, 10)
** 
** datasets:
** biological replicate,filename/URL/{SRR,ERR,DRR} reference
** 
*/

def options = ['bfs':null, 'ltr':null, 'rfs':null, 'dist':null, 'mapq':null, 'intq':null, 'intc':null]

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
   } else if (line =~ /^datasets:/) {
      status = 1
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
   case 1:
     t = line.split(/\s*,\s*/).toList()
     if (t.size() == 2) {
        if (t[1] =~ /^SRR|^ERR|^DRR/) {
           file_getref << tuple([t[-1]],t[0])
        } else {
           log.info "error: dataset iPCR entries must specify 2 files, 2 URL or a GEO reference starting with SRR/ERR/DRR. Error in entry '$line'"
           exit 1
        }
     } else if (t.size() == 3) {
        if (t[1] =~ /^http|^ftp/ && t[2] =~ /^http|^ftp/) {
           file_getref << tuple([t[-2],t[-1]],t[0])
        } else {
           read1 = file("${t[-2]}")
           read2 = file("${t[-1]}")
           if (read1.isFile() && read2.isFile()) {
              datasets << tuple([read1,read2],t[0])
           } else {
              log.info "error: iPCR files not found, in '$line'. (URLs must start with http/ftp)"
              exit 1
           }
        }
     } else {
        log.info "error: incorrect format of iPCR dataset entry '$line' in $params.options"
        exit 1
     }
     break
   }
}

// Parse options
if (!options['bfs'] || !options['ltr'] || !options['dist']) {
   log.info "error: 'bfs', 'ltr' and 'dist' options must be defined in $params.options before 'datasets:'."
   exit 1
}
if (!isInteger(options['dist'])) {
   log.info "error: 'dist' value provided in $params.options must be a valid integer."
}

// Group the same references/URLS to avoid redundant downloads.
file_getref.close()
file_getref.groupTuple().into{gfile_ref}

process getDataset {
   // Process options
   tag "${data[0]}_${data[1]}"
   publishDir path:"${params.out}/datasets/", mode:'symlink'
   // Cluster options
   cpus 1
   memory '2GB'

   input:
   val data from gfile_ref
   output:
   set file('*.fastq.gz'), replicate into datasets
   script:
   replicate = data[1]
   ref = data[0]
   if (ref.size() == 1) {
      if (ref[0] =~ /^SRR|^ERR|^DRR/) {
         """
         fastq-dump --split-files --gzip -A ${ref[0]}
         rm -f ~/ncbi/public/sra/${ref[0]}.sra
         """
      } else {
         log.info "error: only one iPCR dataset specified (must be PE files or sigle GEO reference)!"
      }
   } else if (ref.size() == 2) {
      """
      wget ${ref[0]}
      wget ${ref[1]}
      """
   }
}

/**********************************************************************
**************************** IPCR PIPELINE ****************************
**********************************************************************/


process extractBarcodes {
   // Process options
   tag "${sample}"
   publishDir path:"${params.out}/map/barcodes_raw/", mode:'symlink'
   
   // Custer options
   cpus 1
   memory '4GB'

   input:
   set file(ipcr), sample from datasets
   file script from extract_barcode_src.first()

   output:
   set file("*.bcd"), file("*.fastq"), sample into raw_barcodes
   
   script:
   rfs_arg = options['rfs'] ? "-r ${options['rfs']}" : ""
   py_options = "-b ${options['bfs']} -l ${options['ltr']} ${rfs_arg}"
   """
   python ${script} ${py_options} --barcodes ipcr_${sample}.bcd --reads ipcr_${sample}.fastq ${ipcr[0]} ${ipcr[1]}
   """
}

process clusterBarcodes {
   // Process options
   tag "$sample"
   publishDir path:"${params.out}/map/barcodes_cluster/", mode:'symlink'
   
   // Cluster options
   cpus 8
   memory '16GB'

   input:
   set file(barcodes), file(reads), sample from raw_barcodes
   file script from barcode_seqnames_src.first()
   
   output:
   set file("ipcr_integ_site_*.fastq"), sample into loci_map
   
   script:
   """
   starcode -d ${options['dist']} --seq-id -t ${task.cpus} ${barcodes} > ipcr_${sample}.stc
   python ${script} ipcr_${sample}.stc ${reads} > ipcr_integ_site_${sample}.fastq
   """
}

process mapIntegs {
   // Process options
   tag "$sample"
   publishDir path:"${params.out}/map/mapping/", mode:'symlink'

   // Cluster options
   cpus 8
   memory '16GB'

   input:
   set file(integ_sites), sample from loci_map
   file index_path from bwa_index.first()
   file index_files from index_files.first()

   output:
   set file("*.bam"), sample into raw_integs
   
   script:
   cat_head = integ_sites.getExtension() == 'gz' ? "zcat ${integ_sites}" : "cat ${integ_sites}"

   """
    OPTS_26='-k17 -r1.3 -B2 -O4 -T22 -L0'
    OPTS_36='-k18 -B3 -O5 -T28 -L0'
    OPTS_40='-L0'
    OPTS_50='-L0'
    SEQLEN=\$((\$(${cat_head} | head -2 | tail -1 | wc -c) - 1));
    if [ \$SEQLEN -le 30 ]; then \
      OPTS=\$OPTS_26; \
    elif [ \$SEQLEN -le 39 ]; then \
      OPTS=\$OPTS_36; \
    elif [ \$SEQLEN -le 46 ]; then \
      OPTS=\$OPTS_40; \
    else \
      OPTS=\$OPTS_50; \
    fi;
    echo \$SEQLEN
    echo \$OPTS
    echo ${index_path}
    bwa mem -t ${task.cpus} \$OPTS ${index_path} ${integ_sites} | samtools view -bS - > ipcr_${sample}.bam
   """
}

process filterRecombinants {
   // Process options
   tag "$sample"
   publishDir path:"${params.out}/map/integs", mode:'symlink'

   // Cluster options
   cpus 1
   memory '4GB'
   
   input:
   set file(raw_integs), sample from raw_integs
   file script from best_assignment_src.first()

   output:
   file "*integs*.txt" into integ_list
   file "*recombinant*.txt" into recombinant_list
   
   script:
   py_args = (options['mapq'] ? "-q ${options['mapq']} " : "") + (options['intq'] ? "--min-score ${options['intq']} " : "") + (options['intc'] ? "--min-reads ${options['intc']}" : "")
   """
   python ${script} ${py_args} <(samtools view ${raw_integs}) ipcr_integs_${sample}.txt ipcr_recombinant_${sample}.txt
   """
}

process mergeIntegs {
   // Process options
   tag "${finteg.size()} files"
   publishDir path:"${params.out}/map/", mode:'symlink'

   // Cluster options
   cpus 1
   memory '4GB'
   
   input:
   file(finteg) from integ_list.flatten().toSortedList()

   output:
   file "ipcr_integs.txt" into integ_file
   
   script:
"""
#!/usr/bin/env python
files = ${finteg.collect{'\"'+it+'\"'}}
for fname in files:
  rep = fname.split('_')[2].split('.')[0]
  with open(fname,'r') as fin, open('ipcr_integs.txt','a') as fout:
    for line in fin:
      fout.write(line.rstrip()+'\\t'+rep+'\\n')
"""
}


// Remove duplicate barcodes.
// If duplicates are integrated in the same locus, keep only one copy.
// Except if they come from different replicates.
process removeDuplicates {
   // Process options
   publishDir path: "${params.out}/", mode: 'copy'
   
   // Cluster options
   cpus 1
   memory '4GB'
   
   input:
   file integs from integ_file
   
   output:
   file 'hiv_integrations.txt' into integs_final

   script:
"""
#!/usr/bin/env python
max_dist = 100
integs = {}
with open('${integs}','r') as f:
  for line in f:
    if line[0] == '#' or line[0:4] == 'brcd': continue
    [brcd,chr,locus,strand,reads,mapq,rep] = line.rstrip().split('\\t')
    key = brcd+rep
    if not integs.has_key(key):
      integs[key] = [chr,int(locus),line]
    else:
      # If len is 1, it means it has been already discarded (set to 0).
      if len(integs[key]) == 1: continue
      # Exists but not yet compared.
      else:
        if integs[key][0] == chr and abs(integs[brcd][1]-int(locus)) < max_dist: continue
        else:
          integs[key] = [0]

with open('hiv_integrations.txt','w+') as fout:
  fout.write('brcd\\tchr\\tlocus\\tstrand\\treads\\tmapq\\trep\\n')
  for v in integs.values():
    if len(v) > 1:
      fout.write(v[2])
"""
}
