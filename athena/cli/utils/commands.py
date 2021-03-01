#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019- Ryan L. Collins <rlcollins@g.harvard.edu>
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
@click.option('-x', '--exclusion-list', default=None,
              help='BED file of regions to exclude, based on SV overlap')
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
@click.option('--af-field', 'af_field', default='AF',
              help='INFO field to use for allele frequency filtering. ' +
              '[default: AF]')
@click.option('--keep-infos', 'keep_infos', default=None,
              help='INFO fields to retain in output VCF (comma-separated). ' + 
                   'Will always retain the following: END, CHR2, SVTYPE, SVLEN. ' + 
                   'Specifying "ALL" will keep all INFO. Note: this does _not_ ' + 
                   'filter records based on INFO. [default: Keep no other INFO]')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output with bgzip')
def filtervcf(vcf, out, chroms, xchroms, svtypes, exclusion_list, minAF, maxAF, 
              minAC, maxAC, minAN, filters, minQUAL, maxQUAL, HWE, af_field, 
              keep_infos, bgzip):
    """
    Filter an input VCF
    """
    utils.filter_vcf(vcf, out, chroms, xchroms, svtypes, exclusion_list, 
                     minAF, maxAF, minAC, maxAC, minAN, filters, 
                     minQUAL, maxQUAL, HWE, af_field, keep_infos, bgzip)


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
@click.argument('outfile_all')
@click.option('--include-chroms', 'chroms', default=None, 
              help='Chromosomes to include (comma-separated) ' + 
              '[default: include all chromosomes]')
@click.option('--exclude-chroms', 'xchroms', default=None, 
              help='Chromosomes to exclude (comma-separated) ' + 
              '[default: exclude no chromosomes]')
@click.option('-s', '--step', 'stepsize', default=None, 
              help='Step size of bins. [default: binsize]')
@click.option('-x', '--exclusion-list-all', default=None, multiple=True,
              help='BED file of regions to exclude for all bins, based on bin ' +
              'overlap. This may be specified multiple times.')
@click.option('--exclusion-list-training', 'exclusion_list_train', default=None, 
              multiple=True, help='BED file of regions to exclude for training ' +
              'bins, based on bin overlap. This may be specified multiple times.')
@click.option('--buffer', 'excl_buffer', type=int, default=0,
              help='Pad exclusion list intervals prior to intersection. If multiple ' + 
              'exclusion lists are specified, elements from each exclusion list will be ' + 
              'padded separately.')
@click.option('--exclusion-list-cov', 'excl_cov', default=1e-9, type=float,
              help='Minimum fraction of bin that must be covered by exclusion ' +
              'lists prior to being excluded [default: any overlap].')
@click.option('--training-bins-out', 'outfile_train', default=None, 
              help='Output BED file of bins for mutation rate training ' + 
              '[default: do not output training bins]')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output with bgzip')
def makebins(genome, binsize, outfile_all, outfile_train, stepsize, 
             exclusion_list_all, exclusion_list_train, excl_buffer, excl_cov, 
             chroms, xchroms, bgzip):
    """
    Create sequential bins
    """
    utils.make_bins(genome, binsize, outfile_all, outfile_train, stepsize, 
                    exclusion_list_all, exclusion_list_train, excl_buffer, 
                    excl_cov, chroms, xchroms, bgzip)


# Slice remote genomic data hosted on Google Cloud Storage
@click.command(name='slice-remote')
@click.argument('urls_tsv', type=click.Path(exists=True))
@click.argument('regions_bed', type=click.Path(exists=True))
@click.option('--ref-fasta', 'ref_fasta', help='Path to reference fasta file ' + \
              '(required if any CRAM files are included in urls_tsv)')
@click.option('--suffix', 'local_suffix', default='local_slice', help='Suffix to ' +
              'append to the filenames of local data slices [default: local_slice]')
@click.option('--updated-tsv', 'tsv_out', default='stdout', help='Path to ' +
              'output tsv with updated paths to local slices of remote data ' +
              '[default: stdout]')
def sliceremote(urls_tsv, regions_bed, ref_fasta, local_suffix, tsv_out):
  """
  Localize slices of remote genomic data
  """
  utils.slice_remote(urls_tsv, regions_bed, ref_fasta, local_suffix, tsv_out)


# Plot feature distributions
@click.command(name='feature-hists')
@click.argument('BED', type=click.Path(exists=True))
@click.argument('png_prefix')
@click.option('--ignore-columns', 'skip_cols', type=int, default=3, 
              help='Skip the first N columns for plotting')
@click.option('--transformations-tsv', 'trans_tsv', help='Two-column tsv listing ' + 
              'all transformations to be applied. Will supersede any transformations ' +
              'passed as arguments.')
@click.option('--log-transform', multiple=True, help='List of column names to ' +
              'be log-transformed prior to plotting. Note that the exact ' +
              'transformation is log10(x+max(x/1000)).')
@click.option('--sqrt-transform', multiple=True, help='List of column names to ' +
              'be square root-transformed prior to plotting.')
@click.option('--exp-transform', multiple=True, help='List of column names to ' +
              'be exponential-transformed prior to plotting.')
@click.option('--square-transform', multiple=True, help='List of column names to ' +
              'be square-transformed prior to plotting.')
@click.option('--boxcox-transform', multiple=True, help='List of column names to ' +
              'be Box-Cox power-transformed prior to decomposition. Note that ' + 
              'the exact transformation is performed on x+max(x/1000).')
def featurehists(bed, png_prefix, skip_cols, trans_tsv, log_transform, sqrt_transform,
                 exp_transform, square_transform, boxcox_transform):
  """
  Plot bin annotation distributions
  """

  if trans_tsv is not None:
      trans = utils.dfutils._load_transformations(trans_tsv)
      log_transform = trans.get('log', [])
      sqrt_transform = trans.get('sqrt', [])
      exp_transform = trans.get('exp', [])
      square_transform = trans.get('square', [])
      boxcox_transform = trans.get('boxcox', [])

  utils.feature_hists(bed, png_prefix, skip_cols, log_transform, sqrt_transform,
                      exp_transform, square_transform, boxcox_transform)


# Transform binwise annotations collected with annotate-bins
@click.command(name='transform')
@click.argument('BED_in', type=click.Path(exists=True))
@click.argument('BED_out')
@click.option('--ignore-columns', 'skip_cols', type=int, default=3, 
              help='Skip the first N columns for plotting')
@click.option('--transformations-tsv', 'trans_tsv', help='Two-column tsv listing ' + 
              'all transformations to be applied. Will supersede any transformations ' +
              'passed as arguments.')
@click.option('--log-transform', multiple=True, help='List of column names to ' +
              'be log-transformed prior to plotting. Note that the exact ' +
              'transformation is log10(x+max(x/1000)).')
@click.option('--sqrt-transform', multiple=True, help='List of column names to ' +
              'be square root-transformed prior to plotting.')
@click.option('--exp-transform', multiple=True, help='List of column names to ' +
              'be exponential-transformed prior to plotting.')
@click.option('--square-transform', multiple=True, help='List of column names to ' +
              'be square-transformed prior to plotting.')
@click.option('--boxcox-transform', multiple=True, help='List of column names to ' +
              'be Box-Cox power-transformed prior to decomposition. Note that ' + 
              'the exact transformation is performed on x+max(x/1000).')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output with bgzip')
def transform(bed_in, bed_out, skip_cols, trans_tsv, log_transform, sqrt_transform,
              exp_transform, square_transform, boxcox_transform, bgzip):
  """
  Transform one or more annotations
  """

  if trans_tsv is not None:
      trans = utils.dfutils._load_transformations(trans_tsv)
      log_transform = trans.get('log', [])
      sqrt_transform = trans.get('sqrt', [])
      exp_transform = trans.get('exp', [])
      square_transform = trans.get('square', [])
      boxcox_transform = trans.get('boxcox', [])

  utils.transform_df(bed_in, bed_out, skip_cols, log_transform, sqrt_transform,
                      exp_transform, square_transform, boxcox_transform, bgzip)


# Intersect SVs and bins
@click.command(name='count-sv')
@click.argument('bins', type=click.Path(exists=True))
@click.argument('sv', type=click.Path(exists=True))
@click.argument('outfile')
@click.option('--sv-format', type=click.Choice(['vcf', 'bed']),
              help='File format of SV input. Required.')
@click.option('--comparison', type=click.Choice(['interval', 'breakpoint']),
              help='Specification of SV-to-bin comparison. Required.')
@click.option('--min-bin-coverage', 'min_cov', type=float, default=0.001,
              help='Minimum fraction of bin to be covered by SV before being ' + 
              'counted. Only used with --comparison breakpoint.')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output with bgzip')
def countsv(bins, sv, outfile, sv_format, comparison, min_cov, bgzip):
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

    utils.count_sv(bins, sv, outfile, sv_format, comparison, min_cov, bgzip)


# Pair bins
@click.command(name='pair-bins')
@click.argument('bins', type=click.Path(exists=True))
@click.argument('outfile_all')
@click.option('--max-dist-all', type=int, default=1000000,
              help='Maximum distance to consider during pairing.')
@click.option('--max-dist-training', 'max_dist_train', type=int, default=1000000,
              help='Maximum distance to consider during pairing for training ' +
              'bins only.')
@click.option('-x', '--exclusion-list-all', default=None, multiple=True,
              help='BED file of regions to exclude for all bin pairs, based on ' +
              'pair span overlap. This may be specified multiple times.')
@click.option('--exclusion-list-training', 'exclusion_list_train', default=None, 
              multiple=True, help='BED file of regions to exclude for training ' +
              'bin pairs, based on pair span overlap. This may be specified ' +
              'multiple times.')
@click.option('--buffer', 'excl_buffer', type=int, default=0,
              help='Pad exclusion list intervals prior to intersection. If ' +
              'multiple exclusion lists are specified, elements from each ' +
              'exclusion list will be padded separately.')
@click.option('--training-pairs-out', 'outfile_train', default=None, 
              help='Output BEDPE file of bin pairs for mutation rate training ' + 
              '[default: do not output training bin pairs]')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output with bgzip')
def pairbins(bins, outfile_all, outfile_train, max_dist_all, max_dist_train, 
             exclusion_list_all, exclusion_list_train, excl_buffer, bgzip):
  """
  Create pairs of bins
  """
  utils.pair_bins(bins, outfile_all, outfile_train, max_dist_all, max_dist_train, 
                  exclusion_list_all, exclusion_list_train, excl_buffer, bgzip)

