#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Create bin-pairs for 2D models
"""

import pandas as pd
import pybedtools
from athena.utils.makebins import _buffer_blacklist
from os import path
from athena.utils.misc import bgzip as bgz


def _get_pairs(chrom, start, end, df, maxdist):
    """
    Return a BedTool of bin pairs within maxdist
    """

    hits = df[(df.chrom == chrom) & (df.start >= end) & (df.start <= end + maxdist)]

    orig_bin = '\t'.join([str(x) for x in [chrom, start, end]])

    pairs = ['\t'.join([orig_bin] + [str(i) for i in x]) for x in hits.values]

    return pybedtools.BedTool('\n'.join(pairs), from_string=True)


def _pair_blacklist(bedpe, blacklists, bl_buffer):
    """
    Exclude pairs from BEDPE based on inter-pair span overlapping blacklist
    """

    if isinstance(blacklists, tuple):
        for bl in blacklists:
            xbt = pybedtools.BedTool(bl)
            xbt = xbt.each(_buffer_blacklist, bl_buffer=bl_buffer).merge()
            bedpe = bedpe.pair_to_bed(b=xbt, type='notospan')
    else:
        xbt = pybedtools.BedTool(blacklists)
        xbt = xbt.each(_buffer_blacklist, bl_buffer=bl_buffer).merge()
        bedpe = bedpe.pair_to_bed(b=xbt, type='notospan')

    return bedpe


def pair_bins(bins, outfile_all, outfile_train, max_dist_all, max_dist_train, 
              blacklist_all, blacklist_train, bl_buffer, bgzip):
    """
    Create pairs of bins as BEDPE from input BED
    """

    # Read bins as pd.dataframe
    df = pd.read_csv(bins, sep='\t', usecols=range(3))
    df.columns = ['chrom', 'start', 'end']

    # Create all pairs
    pairs = [_get_pairs(row[0], row[1], row[2], df, max_dist_all) for row in df.values]
    pairs = pairs[0].cat(*pairs[1:], postmerge=False)

    # Filter all pairs
    if blacklist_all is not None:
        pairs = _pair_blacklist(pairs, blacklist_all, bl_buffer)

    # Save pairs to file
    if '.gz' in outfile_all:
        outfile_all = path.splitext(outfile_all)[0]
    pairs.saveas(outfile_all, trackline='\t'.join(['#chrA','startA','endA', 
                                                   'chrB', 'startB', 'endB']))

    if bgzip:
        bgz(outfile_all)

    # Further filter training pairs by BL & dist
    if outfile_train is not None \
    and max_dist_train < max_dist_all:
        pairs = pairs.filter(lambda x: int(x[4]) <= int(x[2]) + max_dist_train)
        if blacklist_train is not None:
            pairs = _pair_blacklist(pairs, blacklist_train, bl_buffer)

    # Save training pairs to file
    if '.gz' in outfile_train:
        outfile_train = path.splitext(outfile_train)[0]
    pairs.saveas(outfile_train, trackline='\t'.join(['#chrA','startA','endA', 
                                                     'chrB', 'startB', 'endB']))

    if bgzip:
        bgz(outfile_train)
