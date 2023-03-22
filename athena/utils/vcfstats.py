#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Gather SV size and spacing distributions from a VCF
"""


import pysam
from sys import stdin
import pybedtools
from numpy import percentile as pct

def vcf_stats(vcf, quantiles):

    # Open connection to input VCF
    if vcf in '- stdin'.split():
        invcf = pysam.VariantFile(sys.stdin) 
    else:
        invcf = pysam.VariantFile(vcf)
    header = invcf.header

    # Collect SV spacing
    vcfbt = pybedtools.BedTool(vcf)
    dists = list(vcfbt.closest(vcfbt, d=True, t='first', N=True).to_dataframe().iloc[:, -1])

    # Print number of records in VCF
    print('\nSV count: {:,}'.format(len(vcfbt)))

    # Print SV spacing quantiles
    print('\nSV spacing quantiles:' +
          '\n---------------------')
    for q in quantiles.split(','):
        qval = pct(dists, 100 * float(q))
        qform = '{}%: {:,} bp'.format(str(100 * float(q)), int(qval))
        print(qform)

    # Collect SV sizes by class
    sizes = {}
    for record in invcf:
        # Sizes
        svtype = record.info['SVTYPE']
        size = record.info['SVLEN']
        if svtype not in sizes.keys():
            sizes[svtype] = [size, ]
        else:
            sizes[svtype].append(size)
    allsizes = sum(sizes.values(), [])

    # Print size quantiles for all SVs
    print('\nSV size quantiles:' +
          '\n------------------')
    for q in quantiles.split(','):
        qval = pct(allsizes, 100 * float(q))
        qform = '{}%: {:,} bp'.format(str(100 * float(q)), int(qval))
        print(qform)
