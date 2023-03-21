#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019- Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Miscellaneous utilities
"""


import subprocess
from numpy import ceil, floor
import os
import gzip
import csv
import pybedtools as pbt
from tempfile import gettempdir
import warnings
from operator import not_
from athena.utils.math import ci_to_sd


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
                'compressed-fasta' : '.fa.gz .fa.gzip .fasta.gz .fasta.gzip'.split(),
                'hic' : ['hic'],
                'gtf' : ['gtf'],
                'compressed-gtf' : ['gtf.gz']}

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


def vcf2bed(vcf, breakpoints=False, add_ci_to_bkpts=True, ci=0.95, z_extend=6):
    """
    Convert a pysam.VariantFile (vcf) to a BED
    """

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
            left_start_sd = ci_to_sd(2 * abs(cipos[0]), ci)
            left_start = int(max([0, floor(start - (z_extend * left_start_sd))]))
            left_end_sd = ci_to_sd(2 * abs(cipos[1]), ci)
            left_end = int(min([ceil(start + 1 + (z_extend * left_end_sd)), chrom_len]))
            first_bp = bp_bed_fmt.format(chrom, left_start, left_end, var_id, 'POS', start, 
                                         round(left_start_sd, 2), round(left_end_sd, 2), z_extend)

            # Format right breakpoint while accounting for CIEND
            chrom_two_len = vcf.header.contigs[chrom_two].length
            if add_ci_to_bkpts and vcf_has_ciend:
                ciend = record.info.get('CIEND', (0, 0))
            else:
                ciend = (0, 0)
            right_start_sd = ci_to_sd(2 * abs(ciend[0]), ci)
            right_start = int(max([0, floor(end - (z_extend * right_start_sd))]))
            right_end_sd = ci_to_sd(2 * abs(ciend[1]), ci)
            right_end = int(min([ceil(end + 1 + (z_extend * left_end_sd)), chrom_len]))
            second_bp = bp_bed_fmt.format(chrom_two, right_start, right_end, var_id, 'END', end, 
                                          round(right_start_sd, 2), round(right_end_sd, 2), z_extend)
            
            new_interval = first_bp + second_bp
        
        else:
            new_interval = '{0}\t{1}\t{2}\t{3}\n'.format(chrom, start, end, var_id)

        intervals = intervals + new_interval

    return pbt.BedTool(intervals, from_string=True)


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

    if os.path.splitext(snv_mus)[1] in '.bgz .gz .gzip'.split():
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
        elif line == '':
            # Have to specifically check for end-of-file when using readline()
            # readline() returns empy string if EOF reached, causing infinite loop here
            break
        else:
            starts.add(int(line.split('\t')[1]))

    bfile.close()

    starts = list(starts)

    if len(starts) < 2:
        msg = 'Not enough entries in {} to automatically determine binsize. Exiting.'
        exit(msg.format(bed_path))

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


def check_header_compliance(filename, query_bt, header_compliance='loose'):
    """
    Ensure all contigs from query_bt are present in header of filename
    If header_compliance = "loose", will search for approximate matches in header

    Returns a path to a BED file with appropriately-formatted query regions
    """

    query_bed = query_bt.fn
    basename = os.path.basename(filename)

    # Make list of contigs for header compliance check
    contigs = set([x.chrom for x in query_bt])

    # Get header
    ftype = determine_filetype(filename)
    if ftype == 'cram':
        get_header_cmd = 'samtools view -H {} > {}.header'
    elif ftype == 'bam':
        get_header_cmd = 'samtools view -H {} > {}.header'
    elif ftype == 'compressed-vcf':
        get_header_cmd = 'tabix -H {} > {}.header'
    else:
        warn = 'Unable to automatically determine header type for {}; ' + \
               'skipping header compliance check.'
        warnings.warn(warn.format(filename))
        return query_bt
    get_header_cmd = get_header_cmd.format(filename, basename)
    ecode = os.system(get_header_cmd)
    if ecode != 0:
        err = 'ERROR: header query returned non-zero exit status for {}'
        exit(err.format(filename))

    # Read list of eligible contigs from header
    with open('{}.header'.format(basename)) as hfile:
        header = [l.rstrip() for l in hfile.readlines()]
    if ftype in 'cram bam'.split():
        header_contigs = [l.split('\t')[1].split(':')[1] for l in header \
                          if l.startswith('@SQ')]
    else:
        header_contigs = [l.split('<ID=')[1].split(',')[0] for l in header \
                          if l.startswith('##contig')]

    # Build map of contig matches
    remap_contigs = {}
    for contig in contigs:
        if contig not in header_contigs:
            if header_compliance == 'strict':
                err = 'ERROR: contig \'{}\' in query file {} not found in ' + \
                      'header of {}'
                exit(err.format(contig, query_bed, filename))
            else:
                if contig.startswith('chr'):
                    if re.sub('^chr', '', contig) in header_contigs:
                        match = re.sub('^chr', '', contig)
                        remap_contigs[contig] = match
                    else:
                        match = None
                else:
                    if 'chr{}'.format(contig) in header_contigs:
                        match = 'chr{}'.format(contig)
                        remap_contigs[contig] = match
                    else:
                        match = None
                if match is None:
                    msg = 'WARNING: contig \'{}\' from query file {} not found ' + \
                          'in header of {}; continuing anyway due to ' + \
                          '--header-compliance \'loose\'.'
                else:
                    msg = 'WARNING: contig \'{}\' from query file {} not found ' + \
                          'in header of {}; matching to \'{}\' instead.'
                    warnings.warn(msg.format(contig, query_bed, filename, match))

    # Remap contigs in regions_bed based on matches in header
    def __remap(interval, remap_contigs):
        ochrom = interval.chrom
        nchrom = remap_contigs.get(ochrom, ochrom)
        if nchrom is not None:
            interval.chrom = nchrom
        return interval
    query_bt = query_bt.each(__remap, remap_contigs=remap_contigs).\
                        saveas('{}.query.bed'.format(basename))

    return query_bt


def header_compliance_cleanup(filename):
    """
    Cleanup temporary files from header compliance check, if any
    """

    basename = os.path.basename(filename)
    for suffix in 'header query.bed'.split():
        fname = '{}.{}'.format(basename, suffix)
        if os.path.exists(fname):
            os.remove(fname)


def check_contig_naming_scheme(bed):
    """
    Check the contig naming scheme of a BED file (passed as string of path)
    """

    # Check first 10 lines of BED for snapshot of contig names
    if 'compressed' in determine_filetype(bed):
        f_bed = gzip.open(bed, 'rt')
    else:
        f_bed = open(bed)
    chroms = set()
    for i in range(10):
        chrom = f_bed.readline().rstrip().split('\t')[0]
        if '#' in chrom:
            continue
        chroms.add(chrom)
    f_bed.close()

    # Infer contig naming scheme
    has_chrom = ['chr' in k for k in chroms]
    if all(has_chrom):
        return 'has_chr'
    elif all(map(not_, has_chrom)):
        return 'no_chr'
    else:
        err = 'ERROR: inconsistent contig naming scheme detected in {}'
        from sys import exit
        exit(err.format(bed))

