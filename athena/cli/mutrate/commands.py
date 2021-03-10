#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""


import click
from athena import mutrate
from athena.utils.misc import bgzip as bgz
from athena.utils.misc import determine_filetype, make_default_bed_header
from os import path
from gzip import GzipFile
from datetime import datetime
import athena.utils.dfutils as dfutils


@click.command(name='annotate-bins')
@click.argument('bins', type=click.Path(exists=True))
@click.argument('outfile')
@click.option('--include-chroms', default=None, 
              help='Chromosomes to include (comma-separated) ' + 
              '[default: include all chromosomes]')
@click.option('-R', '--ranges', default=None,
              help='BED file containing range(s) for bin restriction.')
@click.option('-t', '--track', default=None, multiple=True,
              help='Path to local annotation track to apply to bins. ' +
              'Also supports remote-hosted BigWig URLs.')
@click.option('-u', '--ucsc-track', default=None, multiple=True, 
              help='UCSC table name to annotate. Requires specifying --ucsc-ref. ' +
              'To extract a specific column for --action map-*, append the column ' +
              'name to this argument with a colon (e.g., "recombRate:decodeAvg"). ' +
              'Default columns will be ignored if multiple columns are specified. ' +
              'Separate multiple columns with commas. Conditional column ' +
              'requirements can be specified by adding a comparison symbol to ' +
              'the column entry (e.g., rmsk:repClass=LINE).')
@click.option('-a', '--actions', default=None, help='Action to apply to each ' + 
              'annotation track. Will be applied sequentially to each entry to ' +
              '--track and --ucsc-track, in that order (all --track entries before ' +
              '--ucsc-track entries). Must be specified the same number of times ' + 
              'as --tracks and --ucsc-tracks combined.',
              multiple=True, 
              type=click.Choice(['count', 'count-unique', 'coverage', 'any-overlap', 
                                 'map-mean', 'map-sum', 'map-min', 'map-max']))
@click.option('-n', '--track-names', default=None, help='Column names to assign to ' +
              'each new column in the header of the annotated bins file. ' +
              'Follows the same rules for ordering as --actions.',
              multiple=True)
@click.option('--track-list', default=None, 
              help='List of local tracks to annotate. Must be specified as three-' +
              'column, tab-delimited text file. One row per track. Columns ' +
              'correspond to argumenmts passed as -t, -a, and -n, respectively. ' +
              'Will be added to other tracks directly specified with -t/-a/-n.')
@click.option('--ucsc-list', default=None, 
              help='List of UCSC tracks to annotate. Must be specified as three-' +
              'column, tab-delimited text file. One row per track. Columns ' +
              'correspond to argumenmts passed as -u, -a, and -n, respectively. ' +
              'Will be added to other tracks directly specified with -u/-a/-n.')
@click.option('-r', '--ucsc-ref', default=None, type=click.Choice(['hg18', 'hg19', 'hg38']),
              help='UCSC reference genome to use with --ucsc-tracks.')
@click.option('--fasta', default=None, help='Reference genome fasta file. If ' +
              'supplied, will annotate all bins with nucleotide content. Will ' +
              'also generate fasta index if not already available locally.')
@click.option('--snv-mutrate', 'snv_mus', default=None, help='tsv file with ' +
              'trinucleotide context mutation rates for SNVs. Will be used to ' +
              'annotate SNV mutation rate per bin if --fasta is also provided.')
@click.option('--no-ucsc-chromsplit', is_flag=True, default=False, 
              help='Disable serial per-chromosome queries for UCSC tracks. May ' + 
              'improve annotation speed. Not recommended unless input bin file ' +
              'is small.')
@click.option('--maxfloat', type=int, default=8, 
              help='Maximum precision of floating-point values. [default: 8]')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output with bgzip.')
@click.option('-q', '--quiet', is_flag=True, default=False, 
              help='Silence progress messages.')
def annotatebins(bins, outfile, include_chroms, ranges, track, ucsc_track, actions, 
                 track_names, track_list, ucsc_list, ucsc_ref, fasta, snv_mus,
                 no_ucsc_chromsplit, maxfloat, bgzip, quiet):
    """
    Annotate bins
    """

    # Sanitize & format inputs
    ucsc_chromsplit = not no_ucsc_chromsplit
    track = list(track)
    ucsc_track = list(ucsc_track)
    actions = tuple([a.lower() for a in actions])

    # Parse file with lists of tracks (if provided) and add to track lists
    if track_list is not None:
      supp_tracks, supp_actions, supp_names = mutrate.parse_track_file(track_list)
      track = track + supp_tracks
      n_ucsc_tracks = len(ucsc_track)
      if n_ucsc_tracks > 0:
        actions = tuple(list(actions[:n_ucsc_tracks]) + supp_actions \
                        + list(actions[n_ucsc_tracks:]))
        track_names = tuple(list(track_names[:n_ucsc_tracks]) + supp_names \
                        + list(track_names[n_ucsc_tracks:]))
      else:
        actions = tuple(list(actions) + supp_actions)
        track_names = tuple(list(track_names) + supp_names)

    # Parse file with list of UCSC tracks (if provided and add to track lists)
    if ucsc_list is not None:
      supp_ucsc_tracks, supp_ucsc_actions, supp_ucsc_names = mutrate.parse_track_file(ucsc_list)
      ucsc_track = ucsc_track + supp_ucsc_tracks
      actions = tuple(list(actions) + supp_ucsc_actions)
      track_names = tuple(list(track_names) + supp_ucsc_names)

    # Handle header reformatting
    n_tracks = len(track) + len(ucsc_track)
    if n_tracks != len(track_names):
        from sys import exit
        err = 'INPUT ERROR: Number of supplied track names ({0}) does not ' + \
              'match number of tracks ({1}).'
        exit(err.format(len(track_names), n_tracks))
    if 'compressed' in determine_filetype(bins):
        header = GzipFile(bins).readline().decode('utf-8').rstrip()
    else:
        header = open(bins, 'r').readline().rstrip()
    if not header.startswith('#'):
      msg = 'INPUT WARNING: '
      status_msg = '[{0}] athena annotate-bins: No header line detected. ' + \
                   'Adding default header.'
      print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S')))
      n_extra_cols = len(header.split('\t')) - 3
      header = make_default_bed_header(n_extra_cols)
    newheader = header + '\t' + '\t'.join(list(track_names))
    if fasta is not None:
        newheader = '\t'.join([newheader, 'pct_gc'])
        if snv_mus is not None:
          newheader = '\t'.join([newheader, 'snv_mu'])
    
    # Annotate bins
    newbins = mutrate.annotate_bins(bins, include_chroms, ranges, track, 
                                    ucsc_track, ucsc_ref, actions, fasta,
                                    snv_mus, maxfloat, ucsc_chromsplit, quiet)

    # Save annotated bins
    if 'compressed' in determine_filetype(outfile):
        outfile = path.splitext(outfile)[0]
    newbins.saveas(outfile, trackline=newheader)

    # Bgzip bins, if optioned
    if bgzip:
        bgz(outfile)


@click.command(name='annotate-pairs')
@click.argument('pairs', type=click.Path(exists=True))
@click.argument('outfile')
@click.option('--include-chroms', 'chroms', default=None, 
              help='Chromosomes to include (comma-separated) ' + 
              '[default: include all chromosomes]')
@click.option('-R', '--ranges', default=None,
              help='BED file containing range(s) for pair restriction.')
@click.option('-t', '--track', default=None, multiple=True,
              help='Path to local annotation track to apply to pairs.')
@click.option('-u', '--ucsc-track', default=None, multiple=True, 
              help='UCSC table name to annotate. Requires specifying --ucsc-ref. ' +
              'Default columns will be ignored if multiple columns are specified. ' +
              'Separate multiple columns with commas. Conditional column ' +
              'requirements can be specified by adding a comparison symbol to ' +
              'the column entry (e.g., rmsk:repClass=LINE).')
@click.option('-a', '--actions', default=None, help='Action to apply to each ' + 
              'annotation track. Will be applied sequentially to each entry to ' +
              '--track and --ucsc-track, in that order (all --track entries before ' +
              '--ucsc-track entries). Must be specified the same number of times ' + 
              'as --tracks and --ucsc-tracks combined.',
              multiple=True, 
              type=click.Choice(['count-pairs', 'pairwise-coverage', 'any-pairwise-overlap']))
@click.option('-n', '--track-names', default=None, help='Column names to assign to ' +
              'each new column in the header of the annotated bins file. ' +
              'Follows the same rules for ordering as --actions.',
              multiple=True)
@click.option('--track-list', default=None, 
              help='List of local tracks to annotate. Must be specified as three-' +
              'column, tab-delimited text file. One row per track. Columns ' +
              'correspond to argumenmts passed as -t, -a, and -n, respectively. ' +
              'Will be added to other tracks directly specified with -t/-a/-n.')
@click.option('--ucsc-list', default=None, 
              help='List of UCSC tracks to annotate. Must be specified as three-' +
              'column, tab-delimited text file. One row per track. Columns ' +
              'correspond to argumenmts passed as -u, -a, and -n, respectively. ' +
              'Will be added to other tracks directly specified with -u/-a/-n.')
@click.option('-r', '--ucsc-ref', default=None, type=click.Choice(['hg18', 'hg19', 'hg38']),
              help='UCSC reference genome to use with --ucsc-tracks.')
@click.option('--fasta', default=None, help='Reference genome fasta file. If ' +
              'supplied, will annotate all bins with nucleotide content. Will ' +
              'also generate fasta index if not already available locally.')
@click.option('--binsize', type=int, default=None, help='Size of bins. [default: ' +
              'infer from spacing of coordinates of pairs]')
@click.option('--homology-cutoff', 'homology_cutoffs', type=float, multiple=True, 
              help='Custom cutoffs for minimum pairwise sequence identity when ' +
              'calculating homology-based features. Requires --fasta. May be specified ' +
              'multiple times. [default: 1.0, 0.9]')
@click.option('--no-ucsc-chromsplit', is_flag=True, default=False, 
              help='Disable serial per-chromosome queries for UCSC tracks. May ' + 
              'improve annotation speed. Not recommended unless input bin file ' +
              'is small.')
@click.option('--maxfloat', type=int, default=8, 
              help='Maximum precision of floating-point values. [default: 8]')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output with bgzip.')
@click.option('-q', '--quiet', is_flag=True, default=False, 
              help='Silence progress messages.')
def annotatepairs(pairs, outfile, chroms, ranges, track, ucsc_track, actions, 
                  track_names, track_list, ucsc_list, ucsc_ref, fasta, binsize, 
                  homology_cutoffs, no_ucsc_chromsplit, maxfloat, bgzip, quiet):
    """
    Annotate pairs
    """

    # Sanitize & format inputs
    ucsc_chromsplit = not no_ucsc_chromsplit
    tracks = list(track)
    ucsc_tracks = list(ucsc_track)
    actions = tuple([a.lower() for a in actions])
    if len(homology_cutoffs) > 0:
      homology_cutoffs = list(homology_cutoffs)
    else:
      homology_cutoffs = [1.0, 0.9]

    # Parse file with lists of tracks (if provided) and add to track lists
    if track_list is not None:
      supp_tracks, supp_actions, supp_names = mutrate.parse_track_file(track_list)
      tracks = tracks + supp_tracks
      n_ucsc_tracks = len(ucsc_tracks)
      if n_ucsc_tracks > 0:
        actions = tuple(list(actions[:n_ucsc_tracks]) + supp_actions \
                        + list(actions[n_ucsc_tracks:]))
        track_names = tuple(list(track_names[:n_ucsc_tracks]) + supp_names \
                        + list(track_names[n_ucsc_tracks:]))
      else:
        actions = tuple(list(actions) + supp_actions)
        track_names = tuple(list(track_names) + supp_names)

    # Parse file with list of UCSC tracks (if provided and add to track lists)
    if ucsc_list is not None:
      supp_ucsc_tracks, supp_ucsc_actions, supp_ucsc_names = mutrate.parse_track_file(ucsc_list)
      ucsc_tracks = ucsc_tracks + supp_ucsc_tracks
      actions = tuple(list(actions) + supp_ucsc_actions)
      track_names = tuple(list(track_names) + supp_ucsc_names)

    # Handle header reformatting
    if 'compressed' in determine_filetype(pairs):
        header = GzipFile(pairs).readline().decode('utf-8').rstrip()
    else:
        header = open(pairs, 'r').readline().rstrip()
    if not header.startswith('#'):
      msg = 'INPUT WARNING: '
      status_msg = '[{0}] athena annotate-pairs: No header line detected. ' + \
                   'Adding default header.'
      print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S')))
      n_extra_cols = len(header.split('\t')) - 3
      header = make_default_bed_header(n_extra_cols)
    if len(track_names) > 0:
      newheader = header + '\t' + '\t'.join(list(track_names))
    else:
      newheader = header
    if fasta is not None:
      for k in homology_cutoffs:
        for direction in 'fwd rev'.split():
          newheader += '\t' + 'longest_{}_kmer_{}pct_identity'.format(direction, int(round(100 * k)))


    # Annotate pairs
    newpairs = mutrate.annotate_pairs(pairs, chroms, ranges, tracks, ucsc_tracks, 
                                      actions, track_names, ucsc_ref, fasta, binsize, 
                                      homology_cutoffs, ucsc_chromsplit, maxfloat, quiet)

    # Save annotated bins
    if 'compressed' in determine_filetype(outfile):
        outfile = path.splitext(outfile)[0]
    newpairs.saveas(outfile, trackline=newheader)

    # Bgzip bins, if optioned
    if bgzip:
        bgz(outfile)


@click.command(name='eigen-bins')
@click.argument('bins', type=click.Path(exists=True))
@click.argument('outfile')
@click.option('-e', '--eigenfeatures', 'components', type=int, default=10,
              help='Number of principal components to return.')
@click.option('-I', '--ICA', 'ica', is_flag=True, default=False,
              help='Perform ICA instead of PCA. In development. [default: PCA]')
@click.option('--min-variance', type=float, default=None,
              help='Optional method for specifying number of components to return. ' + 
              'Specify minimum proportion of variance to explain.')
@click.option('--transformations-tsv', 'trans_tsv', help='Two-column tsv listing ' + 
              'all transformations to be applied. Will supersede any transformations ' +
              'passed as arguments.')
@click.option('--log-transform', multiple=True, help='List of column names to ' +
              'be log-transformed prior to decomposition. Note that the exact ' +
              'transformation is log10(x+max(x/1000)).')
@click.option('--sqrt-transform', multiple=True, help='List of column names to ' +
              'be square root-transformed prior to decomposition.')
@click.option('--exp-transform', multiple=True, help='List of column names to ' +
              'be exponential-transformed prior to decomposition.')
@click.option('--square-transform', multiple=True, help='List of column names to ' +
              'be square-transformed prior to decomposition.')
@click.option('--boxcox-transform', multiple=True, help='List of column names to ' +
              'be Box-Cox power-transformed prior to decomposition. Note that ' + 
              'the exact transformation is applied to x+max(x/1000).')
@click.option('--fill-missing', type=str, default='0',
              help='Behavior for filling missing values. Can specify numeric ' + 
              'value to fill all missing cells, or "mean"/"median" to ' +
              'impute on a per-column basis. [0]')
@click.option('--skip-columns', type=int, default=3,
              help='Skip first N columns of input bins. [3]')
@click.option('--maxfloat', type=int, default=8, 
              help='Maximum precision of floating-point values. [default: 8]')
@click.option('--max-components', 'max_pcs', type=int, default=100, 
              help='Maximum number of components to calculate.')
@click.option('-s', '--stats', default=None,
              help='File to write out Eigenfeature stats.')
@click.option('-p', '--prefix', default='eigenfeature', help='String prefix to ' +
              'use when labeling eigenfeatures.')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output BED with bgzip.')
def annodecomp(bins, outfile, ica, components, min_variance, trans_tsv, log_transform, 
               sqrt_transform, exp_transform, square_transform, boxcox_transform,
               fill_missing, skip_columns, maxfloat, max_pcs, stats, prefix, bgzip):
    """
    Eigendecomposition of annotations
    """

    if trans_tsv is not None:
      trans = dfutils._load_transformations(trans_tsv)
      log_transform = trans.get('log', [])
      sqrt_transform = trans.get('sqrt', [])
      exp_transform = trans.get('exp', [])
      square_transform = trans.get('square', [])
      boxcox_transform = trans.get('boxcox', [])

    mutrate.decompose_bins(bins, outfile, ica, components, min_variance, log_transform, 
                           sqrt_transform, exp_transform, square_transform, 
                           boxcox_transform, fill_missing, skip_columns, maxfloat, 
                           max_pcs, stats, prefix, bgzip)
