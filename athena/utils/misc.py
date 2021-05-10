#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019- Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Miscellaneous utilities
"""


import subprocess
from scipy.stats import chi2, norm
from numpy import ceil, floor
import pybedtools
from os import path
import gzip
import csv
import pybedtools as pbt
from tempfile import gettempdir


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


def bedtool_to_genome_file(bt, contig_size=500000000):
    """
    Creates a synthetic genome file for BEDtools operations based on the
    contigs present in an input pbt.BedTool
    Returns: path to genome BED
    """

    gpath = gettempdir() + '/athena_bins_tmp_genome.tsv'
    gfile = open(gpath, 'w')

    for contig in list(dict.fromkeys([f.chrom for f in bt])):
        gfile.write('{}\t{}\n'.format(contig, contig_size))

    gfile.close()

    return gpath


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


def vcf2bed(vcf, breakpoints=False, add_ci_to_bkpts=True, ci=0.95, z_extend=6):
    """
    Convert a pysam.VariantFile (vcf) to a BED
    """

    def _ci_to_sd(ci_width, ci_pct):
        """
        Infer standard deviation from a confidence interval
        """

        # Note: ci_width measures the full width of the CI (lower to upper bound)
        # Thus, the total number of Z-scores covered in this interval is 2 * norm.ppf of the tail
        zscore_width = abs(2 * norm.ppf((1 - ci_pct) / 2))

        return ci_width / zscore_width

    intervals = ''
    bp_bed_fmt = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'

    vcf_has_cipos = 'CIPOS' in vcf.header.info.keys()
    vcf_has_ciend = 'CIEND' in vcf.header.info.keys()
    
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
            # Format left breakpoint while accounting for CIPOS
            chrom_len = vcf.header.contigs[chrom].length
            if add_ci_to_bkpts and vcf_has_cipos:
                cipos = record.info.get('CIPOS', (0, 0))
            else:
                cipos = (0, 0)
            left_start_sd = _ci_to_sd(2 * abs(cipos[0]), ci)
            left_start = int(max([0, floor(start - (z_extend * left_start_sd))]))
            left_end_sd = _ci_to_sd(2 * abs(cipos[1]), ci)
            left_end = int(min([ceil(start + 1 + (z_extend * left_end_sd)), chrom_len]))
            first_bp = bp_bed_fmt.format(chrom, left_start, left_end, var_id, 'POS', start, 
                                         round(left_start_sd, 2), round(left_end_sd, 2), z_extend)

            # Format right breakpoint while accounting for CIEND
            chrom_two_len = vcf.header.contigs[chrom_two].length
            if add_ci_to_bkpts and vcf_has_ciend:
                ciend = record.info.get('CIEND', (0, 0))
            else:
                ciend = (0, 0)
            right_start_sd = _ci_to_sd(2 * abs(ciend[0]), ci)
            right_start = int(max([0, floor(end - (z_extend * right_start_sd))]))
            right_end_sd = _ci_to_sd(2 * abs(ciend[1]), ci)
            right_end = int(min([ceil(end + 1 + (z_extend * left_end_sd)), chrom_len]))
            second_bp = bp_bed_fmt.format(chrom_two, right_start, right_end, var_id, 'END', end, 
                                          round(right_start_sd, 2), round(right_end_sd, 2), z_extend)
            
            new_interval = first_bp + second_bp
        
        else:
            new_interval = '{0}\t{1}\t{2}\t{3}\n'.format(chrom, start, end, var_id)

        intervals = intervals + new_interval

    return pybedtools.BedTool(intervals, from_string=True)


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

def add_names_to_bed(bedtool):
    """
    Add coordinate-based names to a pbt.BedTool
    """

    bt_str = ''

    newrec_fmt = '{0}\t{1}\t{2}\t{0}_{1}_{2}\n'

    for record in bedtool:
        newrec = newrec_fmt.format(record.chrom, record.start, record.end)
        bt_str += newrec
    
    return pbt.BedTool(bt_str, from_string=True)

