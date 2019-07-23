#!/usr/bin/env python
### From issue42_lakhia; published

# -*- coding: latin-1 -*-
'''Take an BED file and list of Alignments and caluclates RPKM for a minimum score'''

EPILOG = '''
For more details:
        %(prog)s --help
'''

import numpy
import pybedtools
from pybedtools.featurefuncs import normalized_to_length
import pysam
import pandas as pd
import argparse
import csv
import math
import os
from sklearn import preprocessing


def get_args():
    parser = argparse.ArgumentParser(
            description=__doc__, epilog=EPILOG,
            formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('-p', '--peaks',
    	help="The peak file used to calculate RPKM.",
        required = True)
    parser.add_argument('-e','--experiments',
    	help="Comma seperate file of experiment name followd by file location.",
        required = True)
    parser.add_argument('-f','--factor',
    	help="Factor that is being analyzed. Output file header",
        required = True)
    parser.add_argument('--minimum',
    	default = 1,
        type=int,
    	help="The minimum RPKM value that at least 1 of experiements has to have per genomic location. Default is --minimum=1")
    args = parser.parse_args()
    return args


def rpkm(peak_file,aln_file,exp_name,columns):
    ''' Return Pandas Dataframe of RPKM value and location'''
    columns.append(exp_name)
    ## RPKM  =   numReads / (geneLength/1000 * totalNumReads/1,000,000 )
    peak_counts = peak_file.multi_bam_coverage(bams=[aln_file])
    total_counts = reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[2]) for l in pysam.idxstats(aln_file).strip().split('\n')])
    rpkm = peak_counts.each(normalized_to_length, 4, float(math.pow(10,9))/total_counts).saveas("test.bed") ### Changed from 3 to 4
    rpkm_df = rpkm.to_dataframe()
    #os.remove('test.bed')
    rpkm_df.columns = columns
    columns.remove(exp_name)
    return rpkm_df


def rpkm_strand(peak_file,aln_file,exp_name,columns):
    ''' Return Pandas Dataframe of RPKM value and location'''
    columns.append(exp_name)
    ## RPKM  =   numReads / (geneLength/1000 * totalNumReads/1,000,000 )
    peak_counts = peak_file.multi_bam_coverage(bams=[aln_file],s=True)
    total_counts = reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[2]) for l in pysam.idxstats(aln_file,split_lines=True)])
    rpkm = peak_counts.each(normalized_to_length, 6, float(math.pow(10,9))/float(total_counts)).saveas("test.bed")
    rpkm_df = rpkm.to_dataframe()
    #os.remove('test.bed')
    rpkm_df.columns = columns
    columns.remove(exp_name)
    return rpkm_df


def main():
    args = get_args()
    peak_file = pybedtools.BedTool(args.peaks)
    experiment_dict = csv.DictReader(open(args.experiments))

    # Create dataframe from peak file
    rpkm_columns = []
    peak_df = peak_file.to_dataframe()
    columns = list(peak_df.columns.values)
    for exp in experiment_dict:
        experiment = exp['experiment']
        aln_file = exp['file']
        rpkm_columns.append(experiment)
        if 'strand' in columns:
            peak_df = pd.merge(peak_df,rpkm_strand(peak_file,aln_file,experiment,columns))
        else:
            peak_df = pd.merge(peak_df,rpkm(peak_file,aln_file,experiment,columns))


    # Seperate for Enhancers with Names and no Names
    # Filter for alteast Min RPKM in experiment condition
    peak_rpkm_only = peak_df[rpkm_columns]
    scaler = preprocessing.StandardScaler(with_mean=False)
    norm = scaler.fit_transform(peak_rpkm_only.values)
    peak_rpkm_std_robust = pd.DataFrame(data=norm, columns=list(peak_rpkm_only.columns.values), index = peak_rpkm_only.index )
    filtered_rpkm = peak_df.loc[(peak_rpkm_std_robust >= args.minimum).any(axis=1)]

    # Write out RPKM matrix
    filtered_peaks = filtered_rpkm[columns]
    filtered_rpkm.to_csv(args.factor + '_filtered_peaks.tsv', header=True, index=None, sep='\t')
    pybedtools.BedTool.from_dataframe(filtered_peaks).saveas(args.factor + '_filtered_peaks.bed')


if __name__ == '__main__':
    main()

