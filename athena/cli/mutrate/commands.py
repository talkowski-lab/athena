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


@click.command(name='annotate-bins')
@click.argument('bins', type=click.Path(exists=True))
@click.argument('outfile')
@click.option('--chroms', default=None, 
              help='Chromosomes to include (comma-separated) ' + 
              '[default: include all chromosomes]')
@click.option('-R', '--ranges', default=None,
              help='BED file containing range(s) for bin restriction.')
@click.option('-t', '--track', default=None, multiple=True,
              help='Path to local annotation track to apply to bins. ' +
              'Also supports remote-hosted BigWig URLs.')
@click.option('-u', '--ucsc-track', default=None, multiple=True, 
              help='UCSC table name to annotate. Can point to either a table or ' + 
              'a linked remote BigWig file. Requires specifying --ucsc-ref. ' +
              'To extract a specific column for --action map-*, append the column ' +
              'name to this argument with a colon (e.g., "recombRate:decodeAvg").')
@click.option('-a', '--actions', default=None, help='Action to apply to each ' + 
              'annotation track. Will be applied sequentially to each entry to ' +
              '--track and --ucsc-track, in that order (all --track entries before ' +
              '--ucsc-track entries). Must be specified the same number of times ' + 
              'as --tracks and --ucsc-tracks combined.',
              multiple=True, 
              type=click.Choice(['count', 'count-unique', 'coverage',
                                 'map-mean', 'map-sum', 'map-min', 'map-max']))
@click.option('-n', '--track-names', default=None, help='Column names to assign to ' +
              'each new column in the header of the annotated bins file. ' +
              'Follows the same rules for ordering as --actions.',
              multiple=True)
@click.option('-r', '--ucsc-ref', default=None, type=click.Choice(['hg18', 'hg19', 'hg38']),
              help='UCSC reference genome to use with --ucsc-tracks.')
@click.option('--fasta', default=None, help='Reference genome fasta file. If ' +
              'supplied, will annotate all bins with nucleotide content. Will ' +
              'also generate fasta index if not already available locally.')
@click.option('--maxfloat', type=int, default=5, 
              help='Maximum precision of floating-point values. [5]')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output with bgzip.')
@click.option('-q', '--quiet', is_flag=True, default=False, 
              help='Silence progress messages.')
def annotatebins(bins, outfile, chroms, ranges, track, ucsc_track, ucsc_ref, 
                 actions, track_names, fasta, maxfloat, bgzip, quiet):
    """
    Annotate bins
    """

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
        header = open(bins, 'r').readline().decode('utf-8').rstrip()
    newheader = header + '\t' + '\t'.join(list(track_names))
    if fasta is not None:
        newheader = '\t'.join([newheader, 'pct_gc'])
    
    # Annotate bins
    newbins = mutrate.annotate_bins(bins, chroms, ranges, list(track), 
                                    list(ucsc_track), ucsc_ref, actions, fasta,
                                    maxfloat, quiet)

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
              help='Number of principal components to return. [0]')
@click.option('--fill-missing', type=str, default='0',
              help='Behavior for filling missing values. Can specify numeric ' + 
              'value to fill all missing cells, or "mean"/"median" to ' +
              'impute on a per-column basis. [0]')
@click.option('--skip-columns', type=int, default=3,
              help='Skip first N columns of input bins. [3]')
@click.option('--maxfloat', type=int, default=5, 
              help='Maximum precision of floating-point values. [5]')
@click.option('-s', '--stats', default=None,
              help='File to write out Eigenfeature stats.')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output BED with bgzip.')
def annodecomp(bins, outfile, components, fill_missing, skip_columns, 
               maxfloat, stats, bgzip):
    """
    Eigendecomposition of annotations
    """
    
    mutrate.decompose_bins(bins, outfile, components, fill_missing, skip_columns,
                           maxfloat, stats, bgzip)
