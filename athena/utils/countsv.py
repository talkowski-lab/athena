#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Count overlap of bins and SVs
"""


import pybedtools
import pysam
from os import path
import gzip
from athena.utils.misc import bgzip as bgz
from athena.utils.misc import vcf2bed


def count_sv(bins_in, sv_in, outfile, sv_format, comparison):

    # Load bins & retain header
    if path.splitext(bins_in)[1] in '.gz .bgz .gzip .bgzip'.split():
        bins_header = gzip.open(bins_in, 'r').readline().decode('utf-8').rstrip()
    else:
        bins_header = open(bins_in, 'r').readline().decode('utf-8').rstrip()
    bins = pybedtools.BedTool(bins_in)
    bins_ncolumns = bins.field_count(1)

    # Set breakpoints flag
    if comparison == 'breakpoint':
        breakpoints = True
    else:
        breakpoints = False

    # Parse input SV file depending on format
    if sv_format == 'vcf':
        vcf = pysam.VariantFile(sv_in)
        sv = vcf2bed(vcf, breakpoints=breakpoints)
    import pdb; pdb.set_trace()
    