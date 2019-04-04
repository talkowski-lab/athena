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
@click.option('-t', '--track', default=None, multiple=True,
              help='Path to local annotation track to apply to bins.')
@click.option('-u', '--ucsc-track', default=None, multiple=True, 
              help='UCSC track to download & annotate. Also requires ' + 
                   'specifying --ucsc-ref and --actions.')
@click.option('-a', '--actions', default=None, help='Action to apply to each ' + 
              'annotation track. Will be applied sequentially to each entry to ' +
              '--track and --ucsc-track, in that order (all --track entries before ' +
              '--ucsc-track entries). Must be specified the same number of times ' + 
              'as tracks plus ucsc-tracks.',
              multiple=True, 
              type=click.Choice(['count', 'count-unique', 'coverage']))
@click.option('-n', '--track-names', default=None, help='Column names to assign to ' +
              'each new column in the header of the annotated bins file. ' +
              'Follows the same rules for ordering as --actions.',
              multiple=True)
@click.option('-r', '--ucsc-ref', default=None, type=click.Choice(['hg18', 'hg19', 'hg38']),
              help='UCSC reference genome to use with --ucsc-tracks.')
@click.option('-z', '--bgzip', is_flag=True, default=False, 
              help='Compress output with bgzip')
def annotatebins(bins, outfile, track, ucsc_track, ucsc_ref, actions, 
                 track_names, bgzip):
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
    header = GzipFile(bins).readline().decode('utf-8').rstrip()
    newheader = header + '\t' + '\t'.join(list(track_names))

    # Annotate bins
    newbins = mutrate.annotate_bins(bins, list(track), list(ucsc_track), 
                                    ucsc_ref, actions)

    # Save annotated bins
    if '.gz' in outfile:
        outfile = path.splitext(outfile)[0]
    newbins.saveas(outfile, trackline=newheader)

    # Bgzip bins, if optioned
    if bgzip:
        bgz(outfile)


@click.command(name='decomp')
def annodecomp():
    """
    Decompose bin annotations
    """
    click.echo('Perform eigenvector decomposition on annotated bins (in dev.)')
