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
from numpy import nan
from datetime import datetime
from athena.mutrate import ucsc
import pyBigWig
from os import path


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
def add_bigwig_track(bins, track, action):
    
    # Load bigWig track
    bigwig = pyBigWig.open(track)

    operation = action.replace('map-', '')

    def _bw_lookup(interval, bigwig, operation='mean'):
        if 'chr' in list(bigwig.chroms().keys())[0]:
            if 'chr' not in interval.chrom:
                interval.chrom = 'chr' + interval.chrom
        if operation == 'sum':
            val = bigwig.stats(interval.chrom, interval.start, interval.end, 'mean')[0]
            if val is not None:
                val = val * len(interval)
        else:
            val = bigwig.stats(interval.chrom, interval.start, interval.end, operation)[0]
        return val

    values = [_bw_lookup(f, bigwig, operation) for f in bins]

    df = pd.read_csv(bins.fn, sep='\t', header=None)
    df['newvals'] = values
    bins = pybedtools.BedTool.from_dataframe(df)

    return bins


# Map a bed to bins
def add_bedgraph_track(bins, track, action):

    # Assumes column to map is last
    if isinstance(track, pybedtools.BedTool):
        track = track.sort().saveas()
    else:
        track = pybedtools.BedTool(track).sort().saveas()

    map_col = track.field_count(1)

    operation = action.replace('map-', '')

    bins = bins.map(track, c=map_col, o=operation)

    return bins


# Clean up long floats in last column of bins
def float_cleanup(bins, maxfloat, start_idx):
    df = pd.read_csv(bins.fn, sep='\t', header=None)
    df.iloc[:, start_idx:] = df.iloc[:, start_idx:].replace('.',nan).apply(pd.to_numeric).round(maxfloat)
    bins = pybedtools.BedTool.from_dataframe(df)

    return bins


# Wrapper function to add a single local track
def add_local_track(bins, track, action, maxfloat, quiet):
    if quiet is False:
        status_msg = '[{0}] athena annotate-bins: Adding track "{1}" ' + \
                     'with action "{2}"'
        print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                track, action))

    if action in 'count count-unique coverage'.split():
        bins = add_bedtool_track(bins, track, action)

    elif 'map-' in action:
        if path.splitext(track)[1] in '.bw .bigwig .bigWig .BigWig'.split():
            bins = add_bigwig_track(bins, track, action)
        else:
            bins = add_bedgraph_track(bins, track, action)

    return bins


# Wrapper function to add a single ucsc track
def add_ucsc_track(bins, db, track, action, query_regions, maxfloat, ucsc_ref, quiet):

    table, columns, conditions = ucsc.parse_table_arg(track)

    # Collect data from UCSC 
    if ucsc.table_exists(db, table):
        status_msg = '[{0}] athena annotate-bins: Adding UCSC track ' + \
                     '"{1}" with action "{2}"'
        print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                table, action))
        ures = ucsc.query_ucsc(bins, table, columns, conditions, 
                               db, action, query_regions)
    else:
        from sys import exit
        err = 'UCSC ERROR: Could not find table "{0}" for reference "{1}"'
        exit(err.format(table, ucsc_ref))

    # Add track to bins
    if action in 'count count-unique coverage'.split():
        bins = add_bedtool_track(bins, ures, action)

    elif 'map-' in action:
        if isinstance(ures, pybedtools.BedTool):
            bins = add_bedgraph_track(bins, ures, action)
        else:
            bins = add_bigwig_track(bins, ures, action)

    return bins


# Annotate nucleotide content from a reference fasta
def add_nuc_content(bins, fasta, maxfloat, quiet):
    bins.saveas()
    pct_gc = [f[4] for f in bins.cut(range(3)).nucleotide_content(fi=fasta)]
    df = pd.read_csv(bins.fn, sep='\t', header=None)
    df['pct_gc'] = pct_gc
    bins = pybedtools.BedTool.from_dataframe(df)

    return bins


# Master bin annotation function
def annotate_bins(bins, chroms, ranges, tracks, ucsc_tracks, ucsc_ref, 
                  actions, fasta, maxfloat, quiet):

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
        bins = bins.filter(lambda x: x.chrom in chrlist).saveas()
    if ranges is not None:
        bins = bins.intersect(range, wa=True).saveas()


    # Get count of columns in original bins
    n_cols_old = bins.field_count()


    # Annotate bins with all local tracks
    track_counter = 0
    if len(tracks) > 0:
        for track in tracks:
            action = actions[track_counter]
            bins = add_local_track(bins, track, action, maxfloat, quiet)
            track_counter += 1


    # Annotate bins with all UCSC tracks
    if len(ucsc_tracks) > 0:
        if quiet is False:
            status_msg = '[{0}] athena annotate-bins: Connecting to UCSC ' + \
                         'Genome Browser database'
            print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                    fasta))
        db = ucsc.ucsc_connect(ucsc_ref)
        query_regions = ucsc.collapse_query_regions(bins).saveas()

        # Iterate over tracks
        for track in ucsc_tracks:
            action = actions[track_counter]
            bins = add_ucsc_track(bins, db, track, action, query_regions, 
                                  maxfloat, ucsc_ref, quiet)
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

        bins = add_nuc_content(bins, fasta, maxfloat, quiet)


    # Clean up long floats
    bins = float_cleanup(bins, maxfloat, start_idx=n_cols_old)


    return bins
