#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Annotate a binned genome
"""


import csv
import pybedtools as pbt
import pandas as pd
from numpy import nan
from datetime import datetime
from athena.mutrate import ucsc
from athena.utils.misc import chromsort, load_snv_mus, snv_mu_from_seq, \
                              determine_filetype, bedtool_to_genome_file, \
                              check_header_compliance, header_compliance_cleanup
import copy
from athena.utils.nuc import get_seqs_from_bt
from athena.utils.dfutils import float_cleanup
import pyBigWig
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


def add_bedtool_track(bins, track, action, header_compliance='loose'):
    """
    Extract feature values (as list) for all bins vs. a single BedTool (or BAM/CRAM)
    """

    if isinstance(track, str):
        ftype = determine_filetype(track)
    else:
        ftype = None

    # Check for header inconsistencies for indexed tracks
    if ftype in 'bam cram compressed-vcf'.split():
        query_bins = copy.deepcopy(bins)
        query_bins = check_header_compliance(track, query_bins, header_compliance)
    else:
        query_bins = bins
    
    if action == 'count':
        if ftype in 'bam cram'.split():
            values = [int(f[-4]) for f in query_bins.coverage(track, sorted=True)]
        else:
            values = [int(f[-1]) for f in query_bins.intersect(track, c=True, wa=True)]

    elif action == 'count-unique':
        gfile = bedtool_to_genome_file(query_bins)
        chroms = set([f.split('\t')[0] for f in open(gfile).readlines()])
        bedtool = pbt.BedTool(track).filter(lambda f: f.chrom in chroms).sort(g=gfile).merge()
        values = [int(f[-1]) for f in query_bins.intersect(bedtool, c=True, wa=True)]

    elif action == 'coverage':
        values = [float(f[-1]) for f in query_bins.coverage(track)]

    elif action == 'any-overlap':
        values = [min([1, int(f[-1])]) for f \
                  in query_bins.intersect(track, c=True, wa=True)]

    else:
        from sys import exit
        exit('INPUT ERROR: --action {0} not recognized.'.format(action))

    if ftype in 'bam cram compressed-vcf'.split():
        header_compliance_cleanup(track)

    return values


def add_bigwig_track(bins, track, action):
    """
    Extract feature values (as list) for all bins from a bigWig or bigBed file
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
            val = nan
        return val

    values = pd.Series([_bw_lookup(f, bigwig, operation) for f in bins]).astype(float)

    return values.tolist()


def add_bedgraph_track(bins, track, action):
    """
    Map feature values (as list) for all bins from a bed/bedgraph
    """

    # Assumes column to map is last and sorts to match bins
    gfile = bedtool_to_genome_file(bins)
    chroms = set([f.split('\t')[0] for f in open(gfile).readlines()])
    if isinstance(track, pbt.BedTool):
        track = track.filter(lambda f: f.chrom in chroms).sort(g=gfile).saveas()
    else:
        track = pbt.BedTool(track).filter(lambda f: f.chrom in chroms).sort(g=gfile).saveas()

    map_col = track.field_count(1)

    operation = action.replace('map-', '')

    if map_col > 3:
        values = pd.Series([f[-1] for f in bins.map(track, c=map_col, o=operation)])
        values = values.replace({'.' : nan}).astype(float).tolist()
    else:
        values = [nan] * len(bins)

    return values


def add_local_track(bins, track, action, quiet):
    """
    Wrapper function to add a single local track
    """

    ftype = determine_filetype(track)

    if quiet is False:
        status_msg = '[{0}] athena annotate-bins: Adding track "{1}" ' + \
                     'with action "{2}"'
        print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                track, action))

    if action in 'count count-unique coverage any-overlap'.split():
        values = add_bedtool_track(bins, track, action)

    elif 'map-' in action:
        if ftype == 'bigwig':
            values = add_bigwig_track(bins, track, action)
        else:
            values = add_bedgraph_track(bins, track, action)

    return values


def add_ucsc_track(bins, db, track, action, query_regions, ucsc_ref, 
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
            contigs = chromsort(list(set([f.chrom for f in bins.each(ucsc._check_hg_compliance)])))
            sub_ures = [ucsc.subquery_ucsc(bins, table, columns, conditions, db, 
                                      action, query_regions, k) for k in contigs]
            if isinstance(sub_ures[0], pbt.BedTool):
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
        values = add_bedtool_track(bins, ures, action)

    elif 'map-' in action:
        if isinstance(ures, pbt.BedTool):
            values = add_bedgraph_track(bins, ures, action)
        else:
            values = add_bigwig_track(bins, ures, action)

    return values


def add_nuc_content(bins, fasta, maxfloat):
    """
    Extract GC content (as list) for all bins from a reference fasta
    """
    
    if 'compressed' in determine_filetype(fasta):
        fasta = GzipFile(fasta)
    
    pct_gc = [float(f[4]) for f in bins.cut(range(3)).nucleotide_content(fi=fasta)]

    return pct_gc


def add_snv_mu(bins, fasta, snv_mus, maxfloat):
    """
    Extract SNV mutation rate (as list) for all bins from a reference fasta
    """

    # Extend all bins by 1bp at start and end (need trinucleotide context for mu)
    bins.saveas()
    def _increment_bin(feat, dist=1, start=True, end=True):
        if start:
            feat.start = max([0, feat.start - dist])
        if end:
            feat.end = feat.end + dist
        return feat
    buffbins = pbt.BedTool(bins).each(_increment_bin).saveas()

    snv_mu_dict = load_snv_mus(snv_mus)

    if 'compressed' in determine_filetype(fasta):
        fasta = GzipFile(fasta)

    values = []
    for seq in get_seqs_from_bt(buffbins, fasta):
        mu = snv_mu_from_seq(seq.rstrip(), snv_mu_dict)
        values.append(mu)

    return values


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


    # Load bins. Note: must read contents from file due to odd utf-8 decoding 
    # behavior for bgzipped BED files
    ftype = determine_filetype(bins)
    if ftype is None:
        ftype = 'unknown'
    if 'compressed' in ftype:
        bins = ''.join(s.decode('utf-8') for s in GzipFile(bins).readlines())
    else:
        bins = ''.join(open(bins, 'r').readlines())
    firstline = bins.split('\n')[0].split('\t')
    if firstline[0].startswith('#'):
        colnames = firstline
    else:
        colnames = None
    n_cols_old = len(firstline)
    bins = pbt.BedTool(bins, from_string=True)


    # Subset bins to specific chromosomes/ranges, if optioned
    if chroms is not None:
        chrlist = chroms.split(',')
        bins = bins.filter(lambda x: x.chrom in chrlist).saveas()
    if ranges is not None:
        bins = bins.intersect(range, wa=True).saveas()


    # Note: more efficient (and stable) when adding many annotations to hold 
    # pd.DataFrame of bins with annotations in memory and convert entire 
    # pd.DataFrame back to pbt.BedTool after adding all annotations as columns
    # This appears to be due to peculiarities in pyBedTools handling of wide BED files
    bins_bt = bins.cut(range(3)).saveas()
    bins_df = bins.to_dataframe(names=colnames, comment='#')


    # Annotate bins with all local tracks
    track_counter = 0
    if len(tracks) > 0:
        for track in tracks:
            action = actions[track_counter]
            bins_df['newtrack_{}'.format(track_counter)] = \
                add_local_track(bins_bt, track, action, quiet)
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
            # Ping db connection is still active (UCSC may timeout over sequential long queries)
            # If UCSC connection has timed out, reopen new connection
            try:
                db.ping(True)
            except:
                try:
                    db.close()
                except:
                    pass
                db = ucsc.ucsc_connect(ucsc_ref)

            # Submit UCSC query
            action = actions[track_counter]
            bins_df['newtrack_{}'.format(track_counter)] = \
                add_ucsc_track(bins_bt, db, track, action, query_regions, 
                               ucsc_ref, ucsc_chromsplit, quiet)
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

        bins_df['pct_gc'] = add_nuc_content(bins, fasta, maxfloat)

        # Annotate bins with SNV mutation rates, if optioned
        if snv_mus is not None:
            if quiet is False:
                status_msg = '[{0}] athena annotate-bins: Adding SNV mutation ' + \
                             'rates from reference fasta "{1}".'
                print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                        fasta))

            bins_df['snv_mu'] = add_snv_mu(bins, fasta, snv_mus, maxfloat)
    

    # Clean up long floats
    bins_df = float_cleanup(bins_df, maxfloat, start_idx=n_cols_old)


    # Return bins as pbt.BedTool
    return pbt.BedTool.from_dataframe(bins_df)

