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


def _apply_exclusion_list(bins, exclusion_lists, excl_buffer, excl_cov):
    if isinstance(exclusion_lists, tuple):
        for bl in exclusion_lists:
            xbt = pybedtools.BedTool(bl)
            xbt = xbt.each(_buffer_exclusion_list, excl_buffer=excl_buffer).merge()
            bins = bins.intersect(xbt, v=True, f=excl_cov)
    else:
        xbt = pybedtools.BedTool(exclusion_lists)
        xbt = xbt.each(_buffer_exclusion_list, excl_buffer=excl_buffer).merge()
        bins = bins.intersect(xbt, v=True, f=excl_cov)
    return bins


def _buffer_exclusion_list(interval, excl_buffer):
    interval.start = max([0, interval.start - excl_buffer])
    interval.stop = interval.stop + excl_buffer
    return interval


def make_bins(genome, binsize, outfile_all, outfile_train, stepsize, 
              exclusion_list_all, exclusion_list_train, excl_buffer, excl_cov, 
              chroms, xchroms, bgzip):

    # Create bins
    if stepsize is None:
        stepsize = binsize
    bins = pybedtools.BedTool().window_maker(g=genome, w=binsize, s=stepsize)

    # Drop bins produced smaller than desired bin size
    # These are sometimes generated by bedtools at q-terminal ends of chromosomes
    bins = bins.filter(lambda x: len(x) >= binsize)

    # Include and/or exclude specified chromosomes
    if chroms is not None:
        bins = bins.filter(lambda f: f.chrom in chroms.split(','))
    if xchroms is not None:
        bins = bins.filter(lambda f: f.chrom not in xchroms.split(','))

    # Filter all bins
    if exclusion_list_all is not None:
        bins = _apply_exclusion_list(bins, exclusion_list_all, excl_buffer, excl_cov)

    # Save bins
    if '.gz' in outfile_all:
        outfile_all = path.splitext(outfile_all)[0]
    bins.saveas(outfile_all, trackline='\t'.join(['#chr','start','end']))

    if bgzip:
        bgz(outfile_all)

    # Also generate training bins, if optioned
    if outfile_train is not None:
        # Create & exclusion_list training bins
        if exclusion_list_train is not None:
            tbins = _apply_exclusion_list(bins, exclusion_list_train, excl_buffer, excl_cov)
        else:
            tbins = bins

        # Save training bins
        if '.gz' in outfile_train:
            outfile_train = path.splitext(outfile_train)[0]
        tbins.saveas(outfile_train, trackline='\t'.join(['#chr','start','end']))
        
        if bgzip:
            bgz(outfile_train)
