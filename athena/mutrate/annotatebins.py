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
from datetime import datetime
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


# Annotate bins from a bigWig or bigBed file
# def add_bigwig_track(bins, track, action):



# Annotate nucleotide content from a reference fasta
def add_nuc_content(bins, fasta, quiet=False):

    bins.saveas()
    pct_gc = [f[4] for f in bins.cut(range(3)).nucleotide_content(fi=fasta)]
    df = pd.read_csv(bins.fn, sep='\t', header=None)
    df['pct_gc'] = pct_gc
    bins = pybedtools.BedTool.from_dataframe(df)
    
    return bins


# Master bin annotation function
def annotate_bins(bins, chroms, ranges, tracks, ucsc_tracks, ucsc_ref, 
                  actions, fasta, quiet):

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


    # Load bins & subset to specific chromosomes/ranges, if optioned
    bins = pybedtools.BedTool(bins)
    if chroms is not None:
        chrlist = chroms.split(',')
        bins = bins.filter(lambda x: x.chrom in chrlist)
    if ranges is not None:
        bins = bins.intersect(range, wa=True)


    # Annotate bins with all local tracks
    track_counter = 0
    if len(tracks) > 0:
        for track in tracks:

            action = actions[track_counter]

            if quiet is False:
                status_msg = '[{0}] athena annotate-bins: Adding track "{1}" ' + \
                             'with action "{2}"'
                print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                        track, action))

            bins = add_bedtool_track(bins, track, action)

            track_counter += 1


    # Annotate bins with all UCSC tracks
    if len(ucsc_tracks) > 0:

        # Get UCSC query regions
        query_regions = ucsc.collapse_query_regions(bins)

        # Connect to UCSC database
        if quiet is False:
            status_msg = '[{0}] athena annotate-bins: Connecting to UCSC ' + \
                         'Genome Browser database'
            print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                    fasta))
        db = ucsc.ucsc_connect(ucsc_ref)

        # Iterate over tracks
        for track in ucsc_tracks:

            action = actions[track_counter]

            # Collect data from UCSC 
            if ucsc.table_exists(db, track):
                status_msg = '[{0}] athena annotate-bins: Adding UCSC track ' + \
                             '"{1}"" with action "{2}"'
                print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                        track, action))
                ures = ucsc.query_ucsc(bins, track, db, action, query_regions)
            else:
                from sys import exit
                err = 'UCSC ERROR: Could not find table "{0}" for reference "{1}"'
                exit(err.format(track, ucsc_ref))

            # Add track to bins
            if action in 'count count-unique coverage'.split():
                bins = add_bedtool_track(bins, ures, action)
            elif 'map-' in action:
                import pdb; pdb.set_trace()

            track_counter += 1

        # Close UCSC connection
        db.close()


    # Annotate bins with nucleotide content, if optioned
    if fasta is not None:

        if quiet is False:
            status_msg = '[{0}] athena annotate-bins: Adding nucleotide ' + \
                         'content from reference fasta "{1}".'
            print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                    fasta))

        bins = add_nuc_content(bins, fasta, quiet)


    return bins
