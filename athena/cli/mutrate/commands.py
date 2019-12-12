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
from os import path
from gzip import GzipFile
from datetime import datetime


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
@click.option('--no-ucsc-chromsplit', is_flag=True, default=False, 
              help='Disable serial per-chromosome queries for UCSC tracks. May ' + 
              'improve annotation speed. Not recommended unless input bin file ' +
              'is small.')
@click.option('--maxfloat', type=int, default=5, 
              help='Maximum precision of floating-point values. [5]')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output with bgzip.')
@click.option('-q', '--quiet', is_flag=True, default=False, 
              help='Silence progress messages.')
def annotatebins(bins, outfile, include_chroms, ranges, track, ucsc_track, ucsc_ref, 
                 actions, track_names, track_list, ucsc_list, fasta, 
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
    if path.splitext(bins)[1] in '.bgz .gz .gzip'.split():
        header = GzipFile(bins).readline().decode('utf-8').rstrip()
    else:
        header = open(bins, 'r').readline().rstrip()
    if not header.startswith('#'):
      msg = 'INPUT WARNING: '
      status_msg = '[{0}] athena annotate-bins: No header line detected. ' + \
                   'Adding default header.'
      print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S')))
      n_addl_col = len(header.split('\t')) - 3
      header = '#chr\tstart\tend'
      if n_addl_col > 0:
        default_colname = 'user_col_{0}'
        default_cols = [default_colname.format(str(i+1)) for i in range(n_addl_col)]
        header = header + '\t' + '\t'.join(default_cols)
    newheader = header + '\t' + '\t'.join(list(track_names))
    if fasta is not None:
        newheader = '\t'.join([newheader, 'pct_gc'])
    
    # Annotate bins
    newbins = mutrate.annotate_bins(bins, include_chroms, ranges, track, 
                                    ucsc_track, ucsc_ref, actions, fasta,
                                    maxfloat, ucsc_chromsplit, quiet)

    # Save annotated bins
    if '.gz' in outfile:
        outfile = path.splitext(outfile)[0]
    newbins.saveas(outfile, trackline=newheader)

    # Bgzip bins, if optioned
    if bgzip:
        bgz(outfile)


@click.command(name='eigen-bins')
@click.argument('bins', type=click.Path(exists=True))
@click.argument('outfile')
@click.option('-e', '--eigenfeatures', 'components', type=int, default=10,
              help='Number of principal components to return. [10]')
@click.option('--min-variance', type=float, default=None,
              help='Optional method for specifying number of PCs to return. ' + 
              'Specify minimum proportion of variance to explain.')
@click.option('--log-transform', multiple=True, help='List of column names to ' +
              'be log-transformed prior to decomposition. Note that the exact ' +
              'transformation is log10(x+max(x/1000)).')
@click.option('--sqrt-transform', multiple=True, help='List of column names to ' +
              'be square root-transformed prior to decomposition.')
@click.option('--fill-missing', type=str, default='0',
              help='Behavior for filling missing values. Can specify numeric ' + 
              'value to fill all missing cells, or "mean"/"median" to ' +
              'impute on a per-column basis. [0]')
@click.option('--skip-columns', type=int, default=3,
              help='Skip first N columns of input bins. [3]')
@click.option('--maxfloat', type=int, default=5, 
              help='Maximum precision of floating-point values. [5]')
@click.option('--max-pcs', type=int, default=100, help='Maximum number of ' +
              'principal components to calculate. [100]')
@click.option('-s', '--stats', default=None,
              help='File to write out Eigenfeature stats.')
@click.option('-p', '--prefix', default='eigenfeature', help='String prefix to ' +
              'use when labeling eigenfeatures.')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output BED with bgzip.')
def annodecomp(bins, outfile, components, min_variance, log_transform, 
               sqrt_transform, fill_missing, skip_columns, maxfloat, max_pcs, 
               stats, prefix, bgzip):
    """
    Eigendecomposition of annotations
    """

    mutrate.decompose_bins(bins, outfile, components, min_variance, log_transform, 
                           sqrt_transform, fill_missing, skip_columns, maxfloat, 
                           max_pcs, stats, prefix, bgzip)
