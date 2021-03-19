#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019- Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Create bin-pairs for 2D models
"""

import pybedtools as pbt
from pysam import TabixFile
from athena.utils import determine_filetype
from athena.utils.exclusionbed import load_exclusion_bts
from os import path
from athena.utils.misc import bgzip as bgz
import gzip


def _get_pairs(outfile, query_vals, bins_tabix, max_dist, xbt, annotate_dist, 
               sort_features, annotate_absdiff, maxfloat):
    """
    Find & curate all candidate pairs for a single input query bin and write to outfile
    """

    # Assign query variables
    chrom = query_vals[0]
    query_start, query_end = [int(k) for k in query_vals[1:3]]
    query_feats = [round(float(k), maxfloat) for k in query_vals[3:]]

    # Use tabix to extract all candidate partners for query bin
    qregion = '{}:{}-{}'.format(chrom, query_start, query_end + max_dist)
    mates = bins_tabix.fetch(qregion)

    # Process each mate in serial
    for mate in mates:
        # Assign mate variables
        mate_vals = mate.split('\t')
        mate_start, mate_end = [int(k) for k in mate_vals[1:3]]
        mate_feats = [round(float(k), maxfloat) for k in mate_vals[3:]]

        # Only keep mates matching or downstream (greater coordinate) of query
        if mate_start < query_start:
            continue

        # Always report query coordinates as (chrom, query_start, mate_end)
        pair_vals = [chrom, query_start, mate_end]

        # Only keep pairs that do not overlap any intervals in xbt
        pair_bt = pbt.BedTool('\t'.join([str(x) for x in pair_vals]), from_string=True)
        if len(pair_bt.intersect(xbt)) > 0:
            continue

        # Compute midpoint distance, if optioned
        if annotate_dist:
            pair_dist = abs(mate_start - query_start)
            pair_vals.append(pair_dist)

        # Order features depending on value of sort_features
        for fl, fr in zip(query_feats, mate_feats):
            if sort_features:
                pair_vals += sorted([fl, fr])
            else:
                pair_vals += [fl, fr]
            # Calculate absolute feature difference, if optioned
            if annotate_absdiff:
                pair_vals.append(round(abs(fl - fr), maxfloat))

        outfile.write('\t'.join([str(x) for x in pair_vals]) + '\n')


def pair_bins(query_bins, all_bins, outfile, max_dist, exclusion_list, excl_buffer, 
              annotate_dist, sort_features, annotate_absdiff, maxfloat, bgzip, 
              input_has_header=True):
    """
    Create pairs of bins from input BED
    """

    # Open connection to infiles & outfile
    if determine_filetype(query_bins) == 'compressed-bed':
        fin = gzip.open(query_bins, 'rt')
    else:
        fin = open(query_bins)
    if input_has_header:
        colnames = [k.replace('#', '') for k in fin.readline().rstrip().split('\t')]
    if all_bins is None:
        bins_tabix = TabixFile(bins)
    else:
        bins_tabix = TabixFile(all_bins)
    xbt = load_exclusion_bts(exclusion_list, excl_buffer)

    # Open connection to output file
    out_ftype, out_ext = determine_filetype(outfile, return_extension=True)
    if 'compressed' in out_ftype:
        outpath = outfile.replace(out_ext, 'bed')
    else:
        outpath = outfile
    fout = open(outpath, 'w')

    # Format header and write to outfile
    header = '#chr start end'.split()
    if annotate_dist:
        header.append('distance')
    for fname in colnames[3:]:
        if sort_features:
            fname_suffixes = ['min', 'max']
        else:
            fname_suffixes = ['left', 'right']
        if annotate_absdiff:
            fname_suffixes.append('absdiff')
        header += ['_'.join([fname, v]) for v in fname_suffixes]
    fout.write('\t'.join(header) + '\n')

    # Identify and curate all pairs for each bin in fin
    for query_line in fin.readlines():
        query_vals = query_line.rstrip().split('\t')
        new_pairs = _get_pairs(fout, query_vals, bins_tabix, max_dist, xbt, 
                               annotate_dist, sort_features, annotate_absdiff,
                               maxfloat)

    # Clean up
    fout.close()
    if bgzip:
        bgz(outpath)

