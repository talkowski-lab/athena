#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Annotate a binned genome
"""


import pybedtools
import pandas as pd
from athena.mutrate import ucsc


# Annotate bins with a single BedTool
def add_bedtool_track(bins, track, action):
    
    if action == 'count':
        bins = bins.intersect(track, c=True, wa=True)

    elif action == 'count-unique':
        bedtool = pybedtools.BedTool(track).sort().merge()
        bins = bins.intersect(bedtool, c=True, wa=True)

    elif action == 'coverage':
        bins.saveas()
        cov = [f[-1] for f in bins.coverage(track)]
        df = pd.read_csv(bins.fn, sep='\t', header=None)
        df['cov'] = cov
        bins = pybedtools.BedTool.from_dataframe(df)

    else:
        from sys import exit
        exit('INPUT ERROR: --action {0} not recognized.'.format(action))

    return bins


# Master bin annotation function
def annotate_bins(bins, tracks, ucsc_tracks, ucsc_ref, actions):

    # Parse & sanity check all track inputs
    n_all_tracks = len(tracks) + len(ucsc_tracks)
    if len(actions) != n_all_tracks:
        from sys import exit
        err = 'INPUT ERROR: Number of actions ({0}) does not match number ' + \
              'of tracks ({1}).'
        exit(err.format(len(actions), n_all_tracks))

    if len(ucsc_tracks) > 0:
        if ucsc_ref is None:
            from sys import exit
            exit('INPUT ERROR: --ucsc-ref must be specified if any UCSC ' + 
                 'tracks are requested.')


    # Load bins
    bins = pybedtools.BedTool(bins)

    # Start counter for tracks
    track_counter = 0

    # Annotate bins with all local tracks
    for track in tracks:
        bins = add_bedtool_track(bins, track, action=actions[track_counter])
        track_counter += 1

    # Annotate bins with all UCSC tracks


    return bins
