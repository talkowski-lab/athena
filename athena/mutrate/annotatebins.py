#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Annotate a binned genome
"""


import csv
import pybedtools
import pandas as pd
from numpy import nan
from datetime import datetime
from athena.mutrate import ucsc
from athena.utils.misc import chromsort, load_snv_mus, snv_mu_from_seq, determine_extension
import pyBigWig
import itertools
from gzip import GzipFile
from os import path


def parse_track_file(infile):
    """
    Parse an input text file with a list of tracks to annotate
    """

    tracks, actions, names = ([], [], [])

    with open(infile) as fin:
        reader = csv.reader(fin, delimiter='\t')
        for track, action, name in reader:
            tracks.append(track)
            actions.append(action.lower())
            names.append(name)

    return tracks, actions, names


def add_bedtool_track(bins, track, action):
    """
    Annotate bins with a single BedTool
    """

    if isinstance(track, str):
        ftype = determine_extension(track)
    else:
        ftype = None
    
    if action == 'count':
        if ftype in 'bam cram'.split():
            bins.saveas()
            cov = [f[-4] for f in bins.coverage(track, sorted=True)]
            df = pd.read_csv(bins.fn, sep='\t', header=None, comment='#')
            df['cov'] = cov
            bins = pybedtools.BedTool.from_dataframe(df)
        else:
            bins = bins.intersect(track, c=True, wa=True)

    elif action == 'count-unique':
        bedtool = pybedtools.BedTool(track).sort().merge()
        bins = bins.intersect(bedtool, c=True, wa=True)

    elif action == 'coverage':
        bins.saveas()
        cov = [f[-1] for f in bins.coverage(track)]
        df = pd.read_csv(bins.fn, sep='\t', header=None, comment='#')
        df['cov'] = cov
        bins = pybedtools.BedTool.from_dataframe(df)

    elif action == 'any-overlap':
        bins.saveas()
        counts = [min([1, int(f[-1])]) for f \
                  in bins.intersect(track, c=True, wa=True)]
        df = pd.read_csv(bins.fn, sep='\t', header=None, comment='#')
        df['overlap'] = counts
        bins = pybedtools.BedTool.from_dataframe(df)

    else:
        from sys import exit
        exit('INPUT ERROR: --action {0} not recognized.'.format(action))

    return bins


def add_bigwig_track(bins, track, action):
    """
    Annotate bins from a bigWig or bigBed file
    """
    
    # Load bigWig track
    bigwig = pyBigWig.open(track)

    operation = action.replace('map-', '')

    def _bw_lookup(interval, bigwig, operation='mean'):
        # print('Starting {0}:{1}-{2}'.format(interval[0], str(interval[1]), str(interval[2])))
        eligible_chroms = list(bigwig.chroms().keys())
        if 'chr' in eligible_chroms[0]:
            if 'chr' not in interval.chrom:
                interval.chrom = 'chr' + interval.chrom
        if interval.chrom in eligible_chroms \
        and interval.end <= bigwig.chroms(interval.chrom):
            if operation == 'sum':
                val = bigwig.stats(interval.chrom, interval.start, interval.end, 'mean')[0]
                if val is not None:
                    val = val * len(interval)
            else:
                val = bigwig.stats(interval.chrom, interval.start, interval.end, operation)[0]
        else:
            val = None
        return val

    values = [_bw_lookup(f, bigwig, operation) for f in bins]

    df = pd.read_csv(bins.fn, sep='\t', header=None, comment='#')
    df['newvals'] = values
    bins = pybedtools.BedTool.from_dataframe(df)

    return bins


def add_bedgraph_track(bins, track, action):
    """
    Map a bed to bins
    """

    # Assumes column to map is last
    if isinstance(track, pybedtools.BedTool):
        track = track.sort().saveas()
    else:
        track = pybedtools.BedTool(track).sort().saveas()

    map_col = track.field_count(1)

    operation = action.replace('map-', '')

    bins = bins.map(track, c=map_col, o=operation)

    return bins


def float_cleanup(bins, maxfloat, start_idx):
    """
    Clean up long floats in last column of bins
    """

    df = pd.read_csv(bins.fn, sep='\t', header=None, comment='#')
    df.iloc[:, start_idx:] = df.iloc[:, start_idx:].replace('.',nan).apply(pd.to_numeric).round(maxfloat)
    bins = pybedtools.BedTool.from_dataframe(df)

    return bins


def add_local_track(bins, track, action, maxfloat, quiet):
    """
    Wrapper function to add a single local track
    """

    ftype = determine_extension(track)

    if quiet is False:
        status_msg = '[{0}] athena annotate-bins: Adding track "{1}" ' + \
                     'with action "{2}"'
        print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                track, action))

    if action in 'count count-unique coverage any-overlap'.split():
        bins = add_bedtool_track(bins, track, action)

    elif 'map-' in action:
        if ftype == 'bigwig':
            bins = add_bigwig_track(bins, track, action)
        else:
            bins = add_bedgraph_track(bins, track, action)

    return bins


def add_ucsc_track(bins, db, track, action, query_regions, maxfloat, ucsc_ref, 
                   chromsplit, quiet):
    """
    Wrapper function to add a single ucsc track
    """

    table, columns, conditions = ucsc.parse_table_arg(track)

    # Collect data from UCSC 
    if ucsc.table_exists(db, table):
        status_msg = '[{0}] athena annotate-bins: Adding UCSC track ' + \
                     '"{1}" with action "{2}"'
        print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                table, action))
        if chromsplit:
            contigs = chromsort(list(set([f.chrom for f in bins])))
            contigs = ['chr' + k for k in contigs]
            sub_ures = [ucsc.subquery_ucsc(bins, table, columns, conditions, db, 
                                      action, query_regions, k) for k in contigs]
            if isinstance(sub_ures[0], pybedtools.BedTool):
                ures = sub_ures[0]
                if len(sub_ures) > 1:
                    ures = ures.cat(*sub_ures[1:], postmerge=False)
            elif isinstance(sub_ures[0], str):
                ures = sub_ures[0]
            else:
                err = 'QUERY RESULTS ERROR: Unknown query return format from ' + \
                      'table "{0}" for reference "{1}"'
                exit(err.format(table, ucsc_ref))
        else:
            ures = ucsc.query_ucsc(bins, table, columns, conditions, 
                                   db, action, query_regions)
    else:
        from sys import exit
        err = 'QUERY ERROR: Could not find table "{0}" for reference "{1}"'
        exit(err.format(table, ucsc_ref))

    # Add track to bins
    if action in 'count count-unique coverage any-overlap'.split():
        bins = add_bedtool_track(bins, ures, action)

    elif 'map-' in action:
        if isinstance(ures, pybedtools.BedTool):
            bins = add_bedgraph_track(bins, ures, action)
        else:
            bins = add_bigwig_track(bins, ures, action)

    return bins


def add_nuc_content(bins, fasta, maxfloat):
    """
    Annotate nucleotide content from a reference fasta
    """
    
    bins.saveas()
    if path.splitext(fasta)[1] in '.bgz .gz .gzip'.split():
        fasta = GzipFile(fasta)
    pct_gc = [f[4] for f in bins.cut(range(3)).nucleotide_content(fi=fasta)]
    df = pd.read_csv(bins.fn, sep='\t', header=None, comment='#')
    df['pct_gc'] = pct_gc
    bins = pybedtools.BedTool.from_dataframe(df)

    return bins


def add_snv_mu(bins, fasta, snv_mus, maxfloat):
    """
    Annotate SNV mutation rate from a reference fasta
    """

    # Extend all bins by 1bp at start and end (need trinucleotide context for mu)
    bins.saveas()
    def _increment_bin(feat, dist=1, start=True, end=True):
        if start:
            feat.start = feat.start - dist
        if end:
            feat.end = feat.end + dist
        return feat
    buffbins = pybedtools.BedTool(bins).each(_increment_bin)

    snv_mu_dict = load_snv_mus(snv_mus)

    if path.splitext(fasta)[1] in '.bgz .gz .gzip'.split():
        fasta = GzipFile(fasta)

    mubins_str = ''
    fseqs = buffbins.sequence(fasta).seqfn
    with open(fseqs) as fin:
        for seqheader, seq in itertools.zip_longest(*[fin]*2):
            coords = seqheader.rstrip().replace('>', '').replace(':', '\t')\
                     .replace('-', '\t').split('\t')
            coords_fmt = '{}\t{}\t{}'.format(coords[0], int(coords[1]) + 1, 
                                             int(coords[2]) - 1)
            mu = snv_mu_from_seq(seq.rstrip(), snv_mu_dict)
            newbin = '\t'.join([coords_fmt, str(mu)])
            mubins_str = '\n'.join([mubins_str, newbin])
    mubins = pybedtools.BedTool(mubins_str, from_string=True)
    
    mudf = mubins.to_dataframe(names=['chrom', 'start', 'end', 'snv_mu'])
    df = pd.read_csv(bins.fn, sep='\t', header=None, comment='#')
    df.columns = ['chrom', 'start', 'end'] + df.columns.tolist()[3:]
    newdf = pd.merge(df, mudf, on=['chrom', 'start', 'end'], how='left', sort=False)

    bins = pybedtools.BedTool.from_dataframe(newdf)

    return bins


def annotate_bins(bins, chroms, ranges, tracks, ucsc_tracks, ucsc_ref, 
                  actions, fasta, snv_mus, maxfloat, ucsc_chromsplit, quiet):
    """
    Master bin annotation function
    """

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
    # Note: must read contents from file due to odd utf-8 decoding behavior for
    # bgzipped BED files
    if path.splitext(bins)[1] in '.bgz .gz .gzip'.split():
        bins = ''.join(s.decode('utf-8') for s in GzipFile(bins).readlines())
    else:
        bins = open(bins, 'r').readlines()
    bins = pybedtools.BedTool(bins, from_string=True)
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
                                  maxfloat, ucsc_ref, ucsc_chromsplit, quiet)
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

        bins = add_nuc_content(bins, fasta, maxfloat)

        # Annotate bins with SNV mutation rates, if optioned
        if snv_mus is not None:
            if quiet is False:
                status_msg = '[{0}] athena annotate-bins: Adding SNV mutation ' + \
                             'rates from reference fasta "{1}".'
                print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                        fasta))

            bins = add_snv_mu(bins, fasta, snv_mus, maxfloat)
    

    # Clean up long floats
    bins = float_cleanup(bins, maxfloat, start_idx=n_cols_old)

    return bins

