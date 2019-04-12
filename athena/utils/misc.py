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


# Bgzip a file
def bgzip(filename):
    subprocess.run(['bgzip', '-f', filename])


# Chi-squared test for deviation from Hardy-Weinberg equilibrium
def hwe_chisq(record):
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


# Convert a pysam.VariantFile (vcf) to a BED
def vcf2bed(vcf, breakpoints=False):

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
