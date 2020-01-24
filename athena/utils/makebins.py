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


def make_bins(genome, binsize, outfile_all, outfile_train, stepsize, 
              blacklist_all, blacklist_train, bl_buffer, xchroms, bgzip):

    # Create bins
    if stepsize is None:
        stepsize = binsize
    bins = pybedtools.BedTool().window_maker(g=genome, w=binsize, s=stepsize)

    # Drop bins produced smaller than desired bin size
    # These are sometimes generated by bedtools at q-terminal ends of chromosomes
    bins = bins.filter(lambda x: len(x) >= binsize)

    # Exclude specified chromosomes
    if xchroms is not None:
        bins = bins.filter(lambda f: f.chrom not in xchroms.split(','))

    # Filter all bins
    if blacklist_all is not None:
        bins = _apply_blacklist(bins, blacklist_all, bl_buffer)

    # Save bins
    if '.gz' in outfile_all:
        outfile_all = path.splitext(outfile_all)[0]
    bins.saveas(outfile_all, trackline='\t'.join(['#chr','start','end']))

    # Bgzip bins, if optioned
    if bgzip:
        bgz(outfile_all)

    # Also generate training bins, if optioned
    if outfile_train is not None:
        # Create & blacklist training bins
        if blacklist_train is not None:
            tbins = _apply_blacklist(bins, blacklist_train, bl_buffer)
        else:
            tbins = bins

        # Save training bins
        if '.gz' in outfile_train:
            outfile_train = path.splitext(outfile_train)[0]
        tbins.saveas(outfile_train, trackline='\t'.join(['#chr','start','end']))

        # Bgzip bins, if optioned
        if bgzip:
            bgz(outfile_train)


def _apply_blacklist(bins, blacklists, bl_buffer):
    if isinstance(blacklists, tuple):
        for bl in blacklists:
            xbt = pybedtools.BedTool(bl)
            xbt = xbt.each(_buffer_blacklist, bl_buffer=bl_buffer).merge()
            bins = bins.intersect(xbt, v=True)
    else:
        xbt = pybedtools.BedTool(blacklists)
        xbt = xbt.each(_buffer_blacklist, bl_buffer=bl_buffer).merge()
        bins = bins.intersect(xbt, v=True)
    return bins


def _buffer_blacklist(interval, bl_buffer):
    interval.start = max([0, interval.start - bl_buffer])
    interval.stop = interval.stop + bl_buffer
    return interval
