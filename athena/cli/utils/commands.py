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


# Adds breakpoint uncertainty annotations to an existing VCF
@click.command(name='breakpoint-confidence')
@click.argument('vcf', type=click.Path(exists=True))
@click.argument('out', type=str)
@click.option('--min-ci', type=int, default=0, help='Width of smallest 95% ' + \
              'confidence interval to be assigned to any variant [default: 0 bp]')
@click.option('--overwrite', is_flag=True, default=False, help='Overwrite existing ' + \
              'CIPOS and CIEND annotations, if any [default: skip records with ' + \
              'existing CIPOS and CIEND annotations]')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output with bgzip')
def breakpointconfidence(vcf, out, min_ci, overwrite, bgzip):
    """
    Annotate breakpoint uncertainty
    """
    print('ATHENA DEV WARNING: breakpoint-confidence is still under development ' + \
          'and currently only supports fixed CI values. This will be updated in a ' + \
          'future version.')
    utils.breakpoint_confidence(vcf, out, min_ci, overwrite, bgzip)


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
@click.option('--header-compliance', 'header_compliance', default='loose', 
              type=click.Choice(['loose', 'strict'], case_sensitive=False),
              help='Behavior when checking consistency between remote file header ' +
              'and regions_bed. If \'loose\', will search for nearest possible ' +
              'contig match (e.g., chr1 vs. 1). If \'strict\', any contigs ' +
              'present in regions_bed not also present in remote file header ' +
              'will cause athena to abort.')
def sliceremote(urls_tsv, regions_bed, ref_fasta, local_suffix, tsv_out, 
                header_compliance):
  """
  Localize slices of remote genomic data
  """
  utils.slice_remote(urls_tsv, regions_bed, ref_fasta, local_suffix, tsv_out, 
                     header_compliance)


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
              'be Box-Cox power-transformed prior to plotting. Note that ' + 
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


# Gather summary statistics from feature distributions
@click.command(name='feature-stats')
@click.argument('BED', type=click.Path(exists=True))
@click.option('-o', '--outfile', default='stdout', 
              help='Output .tsv [default: stdout]')
@click.option('--ignore-columns', 'skip_cols', type=int, default=3, 
              help='Skip the first N columns for plotting')
@click.option('--transformations-tsv', 'trans_tsv', help='Two-column tsv listing ' + 
              'all transformations to be applied. Will supersede any transformations ' +
              'passed as arguments.')
@click.option('--log-transform', multiple=True, help='List of column names to ' +
              'be log-transformed. Note that the exact transformation is ' +
              'log10(x+max(x/1000)).')
@click.option('--sqrt-transform', multiple=True, help='List of column names to ' +
              'be square root-transformed.')
@click.option('--exp-transform', multiple=True, help='List of column names to ' +
              'be exponential-transformed.')
@click.option('--square-transform', multiple=True, help='List of column names to ' +
              'be square-transformed.')
@click.option('--boxcox-transform', multiple=True, help='List of column names to ' +
              'be Box-Cox power-transformed. Note that the exact transformation ' +
              'is performed on x+max(x/1000).')
@click.option('--maxfloat', type=int, default=8, 
              help='Maximum precision of floating-point values. [default: 8]')
def featurestats(bed, outfile, skip_cols, trans_tsv, log_transform, sqrt_transform,
                 exp_transform, square_transform, boxcox_transform, maxfloat):
  """
  Compute feature distributions
  """

  if trans_tsv is not None:
      trans = utils.dfutils._load_transformations(trans_tsv)
      log_transform = trans.get('log', [])
      sqrt_transform = trans.get('sqrt', [])
      exp_transform = trans.get('exp', [])
      square_transform = trans.get('square', [])
      boxcox_transform = trans.get('boxcox', [])

  utils.feature_stats(bed, outfile, skip_cols, log_transform, sqrt_transform,
                      exp_transform, square_transform, boxcox_transform, maxfloat)


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
              'be log-transformed. Note that the exact transformation is ' + 
              'log10(x+max(x/1000)).')
@click.option('--sqrt-transform', multiple=True, help='List of column names to ' +
              'be square root-transformed.')
@click.option('--exp-transform', multiple=True, help='List of column names to ' +
              'be exponential-transformed.')
@click.option('--square-transform', multiple=True, help='List of column names to ' +
              'be square-transformed.')
@click.option('--boxcox-transform', multiple=True, help='List of column names to ' +
              'be Box-Cox power-transformed. Note that the exact transformation ' +
              'is performed on x+max(x/1000).')
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


# Pair bins
@click.command(name='pair-bins')
@click.argument('bins', type=click.Path(exists=True))
@click.argument('outfile')
@click.option('--bin-superset', 'all_bins', help='Master BED file of all bins to ' +
              'be considered as potential pairs. This input is only useful ' +
              'if the BINS positional argument does not contain all bins ' +
              'for a given chromosome.')
@click.option('--max-dist', type=int, default=1000000,
              help='Maximum distance to search for candidate pairs.')
@click.option('-x', '--exclusion-list', default=None, multiple=True,
              help='BED file of regions to exclude based on pair span overlap. ' +
              'This may be specified multiple times.')
@click.option('--excl-buffer', 'excl_buffer', type=int, default=0,
              help='Pad exclusion list intervals prior to intersection. If ' +
              'multiple exclusion lists are specified, elements from each ' +
              'exclusion list will be padded separately.')
@click.option('--annotate-distance', 'annotate_dist', is_flag=True, default=False, 
              help='Add a new feature corresponding to the distance between the ' +
              'midpoints of each bin per pair.')
@click.option('--sort-features', 'sort_features', is_flag=True, default=False, 
              help='Report sorted feature values per pair of bins as (min, max). ' +
              '[default: report as (left, right)]')
@click.option('--annotate-absdiff', 'annotate_absdiff', is_flag=True, default=False, 
              help='Add one new feature per input feature corresponding to the ' +
              'absolute difference per feature between each bin per pair.')
@click.option('--maxfloat', type=int, default=8, 
              help='Maximum precision of floating-point values. [default: 8]')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output with bgzip')
def pairbins(bins, all_bins, outfile, max_dist, exclusion_list, excl_buffer, 
             annotate_dist, sort_features, annotate_absdiff, maxfloat, bgzip):
  """
  Create pairs of bins
  """
  utils.pair_bins(bins, all_bins, outfile, max_dist, exclusion_list, excl_buffer, 
                  annotate_dist, sort_features, annotate_absdiff, maxfloat, bgzip)

