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
@click.command(name='vcf-filter')
@click.argument('vcf', type=click.Path(exists=True))
@click.argument('out', type=str)
@click.option('--include-chroms', 'chroms', default=None, 
              help='Chromosomes to include (comma-separated) ' + 
              '[default: include all chromosomes]')
@click.option('--exclude-chroms', 'xchroms', default=None, 
              help='Chromosomes to exclude (comma-separated) ' + 
              '[default: exclude no chromosomes]')
@click.option('--svtypes', default=None, 
              help='SV classes to include (comma-separated) [default: all SVs]')
@click.option('-x', '--blacklist', default=None,
              help='BED file of regions to exclude, based on SV overlap [default: None]')
@click.option('--minAF', 'minAF', type=float, default=0, 
              help='Minimum allowed allele frequency [default: 0]')
@click.option('--maxAF', 'maxAF', type=float, default=1.0, 
              help='Maximum allowed allele frequency [default: 1.0]')
@click.option('--minAC', 'minAC', type=int, default=0, 
              help='Minimum allowed allele count [default: 0]')
@click.option('--maxAC', 'maxAC', type=int, default=None, 
              help='Maximum allowed allele count [default: None]')
@click.option('--minAN', 'minAN', type=int, default=0, 
              help='Minimum allowed allele number [default: 0]')
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
                   'Specifying "ALL" will keep all INFO. Note: this does _not_ ' + 
                   'filter records based on INFO. [default: Keep no other INFO]')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output with bgzip')
def filtervcf(vcf, out, chroms, xchroms, svtypes, blacklist, minAF, maxAF, 
              minAC, maxAC, minAN, filters, minQUAL, maxQUAL, HWE, 
              keep_infos, bgzip):
    """
    Filter an input VCF
    """
    utils.filter_vcf(vcf, out, chroms, xchroms, svtypes, blacklist, 
                     minAF, maxAF, minAC, maxAC, minAN, filters, 
                     minQUAL, maxQUAL, HWE, keep_infos, bgzip)


# Gathers size distributions per SV class from a VCF
@click.command(name='vcf-stats')
@click.argument('vcf', type=click.Path(exists=True))
@click.option('-q', '--quantiles', default='0.25,0.5,0.75,0.9,0.95,0.99,0.999', 
              help='SV size quantiles to compute (comma-separated floats)')
def vcfstats(vcf, quantiles):
    """
    Get SV size & spacing
    """
    utils.vcf_stats(vcf, quantiles)


# Bin creation
@click.command(name='make-bins')
@click.argument('genome', type=click.Path(exists=True))
@click.argument('binsize', type=int)
@click.argument('outfile')
@click.option('-s', '--step', 'stepsize', default=None, 
              help='Step size of bins. [default: binsize]')
@click.option('-x', '--blacklist', default=None, multiple=True,
              help='BED file of regions to exclude, based on bin overlap. ' + 
              'This may be specified multiple times. [default: None]')
@click.option('--buffer', type=int, default=0,
              help='Pad blacklist intervals prior to intersection. If multiple ' + 
              'blacklists are specified, elements from each blacklist will be ' + 
              'padded separately. [default: None]')
@click.option('--exclude-chroms', 'xchroms', default=None, 
              help='Chromosomes to exclude (comma-separated) ' + 
              '[default: exclude no chromosomes]')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output with bgzip')
def makebins(genome, binsize, outfile, stepsize, blacklist, buffer, xchroms, bgzip):
    """
    Create sequential bins
    """
    utils.make_bins(genome, binsize, outfile, stepsize, blacklist, buffer, xchroms, bgzip)


# Pair bins
@click.command(name='pair-bins')
def pair():
    """
    Pair bins
    """
    click.echo('Pair bins (in dev.)')


# Intersect SVs and bins
@click.command(name='count-sv')
@click.argument('bins', type=click.Path(exists=True))
@click.argument('sv', type=click.Path(exists=True))
@click.argument('outfile')
@click.option('--sv-format', type=click.Choice(['vcf', 'bed']),
              help='File format of SV input. Required.')
@click.option('--comparison', type=click.Choice(['interval', 'breakpoint']),
              help='Specification of SV-to-bin comparison. Required.')
def countsv(bins, sv, outfile, sv_format, comparison):
    """
    Intersect SV and bins
    """

    # Ensure --sv-format is specified
    if sv_format not in 'vcf bed'.split():
      from sys import exit
      if sv_format is None:
        err = 'INPUT ERROR: --sv-format is required. Options: "vcf" or "bed".'
      else:
        err = 'INPUT ERROR: --sv-format "{0}" not recognized. Options: "vcf" or "bed".'
      exit(err.format(sv_format))

    # Ensure --comparison is specified
    if comparison not in 'interval breakpoint'.split():
      from sys import exit
      if comparison is None:
        err = 'INPUT ERROR: --comparison is required. Options: ' + \
              '"interval" or "breakpoint".'
      else:
        err = 'INPUT ERROR: --comparison "{0}" not recognized. Options: ' + \
              '"interval" or "breakpoint".'
      exit(err.format(comparison))

    utils.count_sv(bins, sv, outfile, sv_format, comparison)
