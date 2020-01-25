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


def _get_pairs(chrom, start, end, df, maxdist):
    """
    Return a BedTool of bin pairs within maxdist
    """

    hits = df[(df.chrom == chrom) & (df.start >= start) & (df.start <= end + maxdist)]

    orig_bin = '\t'.join([str(x) for x in [chrom, start, end]])

    pairs = ['\t'.join([orig_bin] + [str(i) for i in x]) for x in hits.values]

    return pybedtools.BedTool('\n'.join(pairs), from_string=True)


def _check_bedpe_span(f_bedpe, bl):
    """
    Check if a single BEDPE-style feature overlaps a blacklist based on pair span
    """

    import pdb; pdb.set_trace()



def _pair_blacklist(bedpe, blacklists, bl_buffer):
    """
    Exclude pairs from BEDPE based on inter-pair span overlapping blacklist
    """

    if isinstance(blacklists, tuple):
        for bl in blacklists:
            xbt = pybedtools.BedTool(bl)
            xbt = xbt.each(_buffer_blacklist, bl_buffer=bl_buffer).merge()
            bedpe = bedpe.filter(_check_bedpe_span, bl=xbt)
    else:
        xbt = pybedtools.BedTool(blacklists)
        xbt = xbt.each(_buffer_blacklist, bl_buffer=bl_buffer).merge()
        bedpe = bedpe.filter(_check_bedpe_span, bl=xbt)
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
        bins = _pair_blacklist(bins, blacklist_all, bl_buffer)


        
    import pdb; pdb.set_trace()


