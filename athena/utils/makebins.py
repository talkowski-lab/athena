#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Segment a reference genome into uniform, sequential bins
"""


import pybedtools
from os import path
from athena.utils.misc import bgzip as bgz


def make_bins(genome, binsize, outfile, stepsize, blacklist, buffer, xchroms, bgzip):

    # Create bins
    bins = pybedtools.BedTool().window_maker(g=genome, w=binsize, s=stepsize)

    # Exclude specified chromosomes
    if xchroms is not None:
        bins = bins.filter(lambda f: f.chrom not in xchroms.split(','))

    # Filter bins
    if blacklist is not None:
        if isinstance(blacklist, tuple):
            for bl in blacklist:
                xbt = pybedtools.BedTool(bl)
                xbt = xbt.each(_buffer_blacklist, buffer=buffer).merge()
                bins = bins.intersect(xbt, v=True)
        else:
            xbt = pybedtools.BedTool(blacklist)
            xbt = xbt.each(_buffer_blacklist, buffer=buffer).merge()
            bins = bins.intersect(xbt, v=True)

    # Save bins
    if '.gz' in outfile:
        outfile = path.splitext(outfile)[0]
    bins.saveas(outfile, trackline='\t'.join(['#chr','start','end']))

    # Bgzip bins, if optioned
    if bgzip:
        bgz(outfile)


def _buffer_blacklist(interval, buffer):
    interval.start = max([0, interval.start - buffer])
    interval.stop = interval.stop + buffer
    return interval