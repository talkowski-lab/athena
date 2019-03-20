#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""


import click
from athena import utils


# VCF filtering
@click.command(name='filter-vcf')
@click.argument('vcf', type=click.Path(exists=True))
@click.argument('out', type=str)
@click.option('--chroms', default=None, 
              help='Chromosomes to include (comma-separated) ' + 
              '[default: all chromosomes]')
@click.option('--svtypes', default=None, 
              help='SV classes to include (comma-separated) [default: all SVs]')
@click.option('--minAF', 'minAF', type=float, default=0, 
              help='Minimum allowed allele frequency [default: 0]')
@click.option('--maxAF', 'maxAF', type=float, default=1.0, 
              help='Maximum allowed allele frequency [default: 1.0]')
@click.option('--minAC', 'minAC', type=int, default=0, 
              help='Minimum allowed allele count [default: 0]')
@click.option('--maxAC', 'maxAC', type=int, default=None, 
              help='Maximum allowed allele count [default: None]')
@click.option('--vcf-filters', 'filters', default='PASS', 
              help='VCF FILTER statuses to be included [default: PASS]')
@click.option('--minQUAL', 'minQUAL', type=int, default=0,
              help='Minimum allowed QUAL score [default: 0]')
@click.option('--maxQUAL', 'maxQUAL', type=int, default=None, 
              help='Maximum allowed QUAL score [default: None]')
@click.option('--pHWE', 'HWE', type=float, default=None, 
              help='Minimum Hardy-Weinberg equilibrium p-value to include ' + 
                   '[default: Do not filter on HWE]')
@click.option('--keep-infos', 'keep_infos', default=None,
              help='INFO fields to retain in output VCF (comma-separated). ' + 
                   'Will always retain the following: END, CHR2, SVTYPE, SVLEN. ' + 
                   'Specifying "ALL" will keep all INFO. [default: Keep no other INFO]')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output with bgzip')
def filtervcf(vcf, out, chroms, svtypes, minAF, maxAF, minAC, maxAC, filters, 
              minQUAL, maxQUAL, HWE, keep_infos, bgzip):
    """
    Filter an input VCF
    """
    utils.filter_vcf(vcf, out, chroms, svtypes, minAF, maxAF, minAC, maxAC, filters, 
                     minQUAL, maxQUAL, HWE, keep_infos, bgzip)


# Bin creation
@click.command(name='make-bins')
def makebins():
    """
    Create bins
    """
    click.echo('Create bins (in dev.)')


# Intersect SVs and bins
@click.command(name='count-sv')
def countsv():
    """
    Intersect SV and bins
    """
    click.echo('Intersect SV and bins (in dev.)')
