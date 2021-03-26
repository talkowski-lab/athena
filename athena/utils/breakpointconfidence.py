#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Annotate a VCF with breakpoint uncertainty
"""


from numpy import ceil
import pysam
from athena.utils import determine_filetype, bgzip as bgz
from sys import stdin, stdout
from os import path


def lookup_cis(record):
    """
    Currently a dummy function that returns all zeroes
    TODO: eventually will look up CIPOS and CIEND values from a pre-trained breakpoint
    uncertainty model
    """

    return (0, 0), (0, 0)


def add_uncertainty(record, min_ci=0, overwrite=False):
    """
    Annotate uncertainty for a single variant
    """

    # TODO: look up pre-trained breakpoint uncertainty model from input .json
    # Eventually this will be stored in the two-element tuples cipos and ciend
    # For now, just assigning all of those variables to zero
    cipos, ciend = lookup_cis(record)

    # Expand CIs if they are smaller than min_ci
    min_ci_half = ceil(min_ci / 2)
    if cipos[1] - cipos[0] < min_ci:
        cipos = (-min_ci_half, min_ci_half)
    if ciend[1] - ciend[0] < min_ci:
        ciend = (-min_ci_half, min_ci_half)

    # Annotate CIs
    record.info['IMPRECISE'] = True
    if overwrite or 'CIPOS' not in record.info.keys():
        record.info['CIPOS'] = cipos
    if overwrite or 'CIEND' not in record.info.keys():
        record.info['CIEND'] = ciend

    return record


def breakpoint_confidence(vcf, out, min_ci, overwrite, bgzip):
    """
    Main function for breakpoint uncertainty annotation
    """

    # Open connection to input VCF
    if vcf in '- stdin'.split():
        invcf = pysam.VariantFile(sys.stdin) 
    else:
        invcf = pysam.VariantFile(vcf)
    header = invcf.header

    # Add CIPOS and CIEND definitions to header
    add_to_header = [
        '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
        '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="95% confidence interval around POS for imprecise variants">',
        '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="95% confidence interval around END for imprecise variants">',
    ]
    for line in add_to_header:
        header.add_line(line)

    # Open connection to output VCF
    if out in '- stdout'.split():
        outvcf = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        if 'compressed' in determine_filetype(out):
            out = path.splitext(out)[0]
        outvcf = pysam.VariantFile(out, 'w', header=header)

    # Process each record
    for record in invcf.fetch():
        newrecord = add_uncertainty(record, min_ci, overwrite)
        outvcf.write(newrecord)
    outvcf.close()

    # Bgzip output VCF, if optioned
    if bgzip:
        bgz(out)

