#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019- Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Miscellaneous utilities
"""


import subprocess
from scipy.stats import chi2
import pybedtools
from os import path
import gzip
import csv
import pybedtools as pbt


def bgzip(filename):
    """
    Bgzip a file
    """

    subprocess.run(['bgzip', '-f', filename])


def determine_filetype(path, return_extension=False):
    """
    Determine file extension for common genomic data formats
    """

    # Enumerate candidate suffix matches
    suf_dict = {'cram' : ['cram'],
                'bam' : ['bam'],
                'vcf' : ['vcf'],
                'compressed-vcf' : 'vcf.gz vcf.bgz vcf.gzip vcf.bgzip'.split(),
                'bed' : ['bed'],
                'compressed-bed' : 'bed.gz bed.bgz bed.gzip bed.bgzip'.split(),
                'bigwig' : '.bw .bigwig .bigWig .BigWig'.split(),
                'fasta' : '.fa .fasta'.split(),
                'compressed-fasta' : '.fa.gz .fa.gzip .fasta.gz .fasta.gzip'.split()}

    for ftype, suffs in suf_dict.items():
        suf_hits = [s for s in suffs if path.endswith(s)]
        if len(suf_hits) > 0:
            if return_extension:
                return ftype, suf_hits[0]
            else:
                return ftype

    # If no matches are found, return None
    if return_extension:
        return None, None
    else:
        return None


def chromsort(contigs):
    """
    Sort a list of strings according to chromosome order
    """

    def _clean_contig_prefixes(contig):
        return contig.replace('chr', '').replace('Chr', 'chr')

    def _is_numeric_contig(contig):
        contig = _clean_contig_prefixes(contig)
        try:
            int(contig)
            return True
        except:
            return False

    def numeric_sort(contigs):
        return sorted(contigs, key=lambda k: int(_clean_contig_prefixes(k)))

    nc = numeric_sort([k for k in contigs if _is_numeric_contig(k)])
    nnc = sorted([k for k in contigs if not _is_numeric_contig(k)])

    return nc + nnc


def hwe_chisq(record):
    """
    Chi-squared test for deviation from Hardy-Weinberg equilibrium
    """

    # Get observed genotype counts
    aa_obs = record.info['N_HOMREF']
    Aa_obs = record.info['N_HET']
    AA_obs = record.info['N_HOMALT']
    total_N = aa_obs + Aa_obs + AA_obs

    # Get expected genotype counts
    a_freq = ((2 * aa_obs) + Aa_obs) / (2 * total_N)
    A_freq = ((2 * AA_obs) + Aa_obs) / (2 * total_N)
    aa_exp = (a_freq ** 2) * total_N
    Aa_exp = 2 * a_freq * A_freq * total_N
    AA_exp = (A_freq ** 2) * total_N

    # Calculate chi-square stats for each genotype
    aa_chisq = ((aa_obs - aa_exp) ** 2) / aa_exp
    Aa_chisq = ((Aa_obs - Aa_exp) ** 2) / Aa_exp
    AA_chisq = ((AA_obs - AA_exp) ** 2) / AA_exp
    chisq = aa_chisq + Aa_chisq + AA_chisq

    p = 1 - chi2.cdf(chisq, 1)

    return p


def vcf2bed(vcf, breakpoints=False):
    """
    Convert a pysam.VariantFile (vcf) to a BED
    """

    intervals = ''
    
    for record in vcf:

        chrom = record.chrom
        if 'CHR2' in record.info.keys():
            chrom_two = record.info['CHR2']
        else:
            chrom_two = chrom
        start = record.pos
        end = record.stop
        var_id = record.id
        
        if breakpoints == True:
            first_bp = '{0}\t{1}\t{2}\t{3}\n'.format(chrom, start, start + 1, var_id)
            second_bp = '{0}\t{1}\t{2}\t{3}\n'.format(chrom_two, end, end + 1, var_id)
            new_interval = first_bp + second_bp
        else:
            new_interval = '{0}\t{1}\t{2}\t{3}\n'.format(chrom, start, end, var_id)

        intervals = intervals + new_interval

    bed = pybedtools.BedTool(intervals, from_string=True)

    return bed


def make_default_bed_header(n_extra_cols):
    """
    Create default BED3+ header line for files lacking informative headers
    """

    header = '#chr\tstart\tend'

    if n_extra_cols > 0:
        default_colname = 'user_col_{0}'
        default_cols = [default_colname.format(str(i+1)) for i in range(n_extra_cols)]
        header = header + '\t' + '\t'.join(default_cols)

    return header


def load_snv_mus(snv_mus):
    """
    Collapse all possible SNV mutation rates per trinucleotide context
    """

    snv_mu_dict = {}

    if path.splitext(snv_mus)[1] in '.bgz .gz .gzip'.split():
        fin = gzip.open(snv_mus, 'rt')
    else:
        fin = open(snv_mus)
    next(fin)
        
    tsv = csv.reader(fin, delimiter='\t')

    for ref, alt, rate in tsv:
        if ref in snv_mu_dict.keys():
            snv_mu_dict[ref] += float(rate)
        else:
            snv_mu_dict[ref] = float(rate)

    return snv_mu_dict


def snv_mu_from_seq(seq, snv_mu_dict):
    """
    Calculate total SNV mutation rate for a given DNA sequence
    """

    mus = [snv_mu_dict.get(seq[(i-1):(i+2)], 0) for i in range(1, len(seq)-1)]

    return sum(mus)


def calc_binsize(bed_path, sample_n_starts=20):
    """
    Estimates bin size from an athena BED based on the minimal distance between the
    start coordinates of the first $sample_starts records with different start 
    positions in a BED file (this assumes the BED file is coordinate-sorted)
    """

    starts = set()

    if determine_filetype(bed_path) == 'compressed-bed':
        bfile = gzip.open(bed_path, 'rt')
    else:
        bfile = open(bed_path)

    while len(starts) < sample_n_starts:
        line = bfile.readline().rstrip()
        if line.startswith('#'):
            continue
        else:
            starts.add(int(line.split('\t')[1]))

    bfile.close()

    starts = list(starts)

    dists = set()
    for i in range(len(starts)):
        for k in range(i+1, len(starts)):
            dists.add(abs(starts[k] - starts[i]))

    return min(list(dists))

