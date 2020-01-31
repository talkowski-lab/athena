#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
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


def bgzip(filename):
    """
    Bgzip a file
    """

    subprocess.run(['bgzip', '-f', filename])


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

