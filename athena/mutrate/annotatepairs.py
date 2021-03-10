#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Annotate pairs of bins
"""

from athena.mutrate import ucsc
from athena.utils.misc import chromsort, calc_binsize, determine_filetype
from athena.utils import nuc
from athena.utils.dfutils import float_cleanup
from gzip import GzipFile
import pybedtools as pbt
from os import path, system
import pysam
from numpy import nan
import pandas as pd
from datetime import datetime


def _split_pairs(pair_interval, binsize):
    """
    Splits a pair (represented as pbt.Interval) into its constituent bins
    """

    chrom, left, right = pair_interval.fields[0:3]
    
    lbin = '\t'.join([chrom, left, str(int(left) + binsize)])
    rbin = '\t'.join([chrom, str(int(right) - binsize), right])

    return pbt.BedTool('\n'.join([lbin, rbin]), from_string=True)


def _pairs_bed_to_bedpe(pairs_bt, binsize):
    """
    Convert athena pairs BED to BEDPE format. Returns pbt.BedTool
    """

    bedpe = ''

    for pair in pairs_bt:
        chrom = pair.chrom
        left_start = pair.start
        left_end = pair.start + binsize
        right_start = pair.end - binsize
        right_end = pair.end
        bedpe += '\t'.join([str(x) for x in [chrom, left_start, left_end,
                                             chrom, right_start, right_end]]) + '\n'

    return pbt.BedTool(bedpe, from_string=True)


def _split_pairtopair_by_binpairs(track_hits, pairs_bedpe_bt, counts_only=False):
    """
    Partition a single pbt.BedTool into a list of pbt.BedTool objects (one per original bin-pair)
    """

    if counts_only:
        starting_value = 0
    else:
        starting_value = ''
    split_res = {'_'.join([pair.chrom, str(pair.start), str(pair[5])]) : starting_value 
                 for pair in pairs_bedpe_bt}

    for hit in track_hits:
        track_name = '_'.join([hit.chrom, str(hit.start), str(hit[5])])
        if counts_only:
            split_res[track_name] += 1
        else:
            split_res[track_name] += '\t'.join(hit[6:12]) + '\n'

    if not counts_only:
        split_res = {tn : pbt.BedTool(h, from_string=True) for tn, h in split_res.items()}

    return split_res


def _bam_to_bedpe(bam_path, query_regions=None, clean_nsorted_copy=True):
    """
    Convert BAM/CRAM records to BEDPE
    """

    # Must namesort BAM/CRAM before running bedtools bamtobed
    sorted_path = path.splitext(bam_path)[0] + '.nsort.bam'
    pysam.sort('-n', '-O', 'bam', '-o', sorted_path, bam_path)

    bedpe = pbt.BedTool(sorted_path).bam_to_bed(bedpe=True).sort().cut(range(6)).saveas()

    # Clean up namesorted BAM/CRAM because it's no longer needed outside of this function
    if clean_nsorted_copy:
        system('rm -rf ' + sorted_path)
    
    return bedpe


def add_pairwise_bedtool_track(pairs_bedpe_bt, track, action, binsize):
    """
    Extract feature values (as list) for all pairs vs. a BedTool (as BEDPE)
    """
    
    # Multiple actions all start with computing bedtools pairtopair -a both
    if action in 'count-pairs pairwise-coverage any-pairwise-overlap'.split():

        track_hits = pairs_bedpe_bt.pairtopair(b=track, type='both')
        if action in 'count-pairs any-pairwise-overlap'.split():
            counts_only = True
        else:
            counts_only = False
        hits_per_bin = _split_pairtopair_by_binpairs(track_hits, pairs_bedpe_bt, 
                                                         counts_only=counts_only)

        if action == 'count-pairs':
            values = list(hits_per_bin.values())

        elif action == 'any-pairwise-overlap':
            values = [min([1, k]) for k in hits_per_bin.values()]

        elif action == 'pairwise-coverage':
            values = []
            for pair, hits in hits_per_bin.items():
                if len(hits) > 0:
                    fields = pair.split('_')
                    pair_interval = pbt.Interval(fields[0], int(fields[1]), int(fields[2]))
                    pairbt = _split_pairs(pair_interval, binsize)
                    covdf_names = 'chr start end items bp total frac'.split()
                    covdf = pairbt.coverage(hits).to_dataframe(names=covdf_names)
                    values.append(covdf.bp.sum() / covdf.total.sum())
                else:
                    values.append(0)

    else:
        from sys import exit
        exit('INPUT ERROR: --action {0} not recognized.'.format(action))

    return values


def add_pairwise_local_track(pairs_bedpe_bt, track, action, query_regions, binsize, quiet):
    """
    Wrapper function to extract values for a single local track
    """

    ftype = determine_filetype(track)

    if quiet is False:
        status_msg = '[{0}] athena annotate-pairs: Adding track "{1}" ' + \
                     'with action "{2}"'
        print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                track, action))

    # Convert BAM/CRAM records to BEDPE, if necessary
    if ftype in 'bam cram'.split():
        track = _bam_to_bedpe(track, query_regions)

    if action in 'count-pairs pairwise-coverage any-pairwise-overlap'.split():
        values = add_pairwise_bedtool_track(pairs_bedpe_bt, track, action, binsize)

    return values


def add_pairwise_ucsc_track(pairs_bedpe_bt, db, track, action, query_regions, 
                            binsize, ucsc_ref, chromsplit, quiet):
    """
    Annotate pairs with a single UCSC track
    """

    table, columns, conditions = ucsc.parse_table_arg(track)

    # Collect data from UCSC 
    if ucsc.table_exists(db, table):
        status_msg = '[{0}] athena annotate-pairs: Adding UCSC track ' + \
                     '"{1}" with action "{2}"'
        print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                table, action))
        if chromsplit:
            contigs = chromsort(list(set([f.chrom for f in pairs_bedpe_bt.each(ucsc._check_hg_compliance)])))
            sub_ures = [ucsc.subquery_ucsc(pairs_bedpe_bt, table, columns, conditions, db, 
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
            ures = ucsc.query_ucsc(pairs_bedpe_bt, table, columns, conditions, 
                                   db, action, query_regions)
    else:
        from sys import exit
        err = 'QUERY ERROR: Could not find table "{0}" for reference "{1}"'
        exit(err.format(table, ucsc_ref))

    # Add track to pairs
    if action in 'count-pairs pairwise-coverage'.split():
        values = add_pairwise_bedtool_track(pairs_bedpe_bt, ures, action, binsize)

    else:
        from sys import exit
        exit('INPUT ERROR: --action {0} not recognized.'.format(action))


    return values


def add_homology(pairs_bt, fasta, binsize, identity=1.0, reverse_seq2=False):
    """
    Compute homology features for all pairs from a reference fasta
    """

    values = []

    # Process pairs in serial
    for pair in pairs_bt:

        # Only keep pairs representing distinct bins
        if pair.length > binsize:

            # Get nucleotide sequences corresponding to both bins in the pair
            bins_bt = _split_pairs(pair, binsize)
            seql, seqr = nuc.get_seqs_from_bt(bins_bt, fasta)

            if reverse_seq2:
                seqr = seqr[::-1]

            longest = nuc.longest_subseq_min_identity(seql, seqr, identity)
            values.append(longest)

        else:
            values.append(nan)

    return values


def annotate_pairs(pairs, chroms, ranges, tracks, ucsc_tracks, actions, track_names, 
                   ucsc_ref, fasta, binsize, homology_cutoffs, ucsc_chromsplit, 
                   maxfloat, quiet):
    """
    Master pair annotation function
    """

    # Infer binsize and filetype
    if binsize is None:
        binsize = calc_binsize(pairs)
    ftype = determine_filetype(pairs)


    # Load pairs. Note: must read contents from file due to odd utf-8 decoding 
    # behavior for bgzipped BED files with pybedtools
    if 'compressed' in ftype:
        pairs = ''.join(s.decode('utf-8') for s in GzipFile(pairs).readlines())
    else:
        pairs = open(pairs, 'r').readlines()
    firstline = pairs.split('\n')[0].split('\t')
    if firstline[0].startswith('#'):
        colnames = firstline
    else:
        colnames = None
    n_cols_old = len(firstline)
    pairs = pbt.BedTool(pairs, from_string=True)


    # Subset pairs to specific chromosomes/ranges, if optioned
    if chroms is not None:
        chrlist = chroms.split(',')
        pairs = pairs.filter(lambda x: x.chrom in chrlist).saveas()
    if ranges is not None:
        pairs = pairs.intersect(range, wa=True).saveas()


    # Note: more efficient (and stable) when adding many annotations to hold 
    # pd.DataFrame of pairs with annotations in memory and convert entire 
    # pd.DataFrame back to pbt.BedTool after adding all annotations as columns
    # This appears to be due to peculiarities in pyBedTools handling of wide BED files
    pairs_bt = pairs.cut(range(3)).saveas()
    pairs_df = pairs.to_dataframe(names=colnames, comment='#')
    pairs_bedpe_bt = _pairs_bed_to_bedpe(pairs_bt, binsize)


    # Make master pbt.BedTool of all bins from all pairs
    split_pair_bts = [_split_pairs(p, binsize) for p in pairs_bt]
    allbins_bt = split_pair_bts[0].cat(*split_pair_bts[1:], postmerge=False).sort().merge(d=-1)
    query_regions = ucsc.collapse_query_regions(allbins_bt).saveas()


    # Annotate bins with all local tracks
    track_counter = 0
    if len(tracks) > 0:
        for track in tracks:
            action = actions[track_counter]
            pairs_df['newtrack_{}'.format(track_counter)] = \
                add_pairwise_local_track(pairs_bedpe_bt, track, action, query_regions, 
                                         binsize, quiet)
            track_counter += 1


    # Annotate bins with all UCSC tracks
    if len(ucsc_tracks) > 0:
        if quiet is False:
            status_msg = '[{0}] athena annotate-pairs: Connecting to UCSC ' + \
                         'Genome Browser database'
            print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                    fasta))
        db = ucsc.ucsc_connect(ucsc_ref)

        # Iterate over tracks
        for track in ucsc_tracks:
            action = actions[track_counter]
            pairs_df['newtrack_{}'.format(track_counter)] = \
                add_pairwise_ucsc_track(pairs_bedpe_bt, db, track, action, query_regions, 
                                        binsize, ucsc_ref, ucsc_chromsplit, quiet)
            track_counter += 1

        # Close UCSC connection
        db.close()


    # Annotate pairs based on nucleotide content, if optioned
    if fasta is not None:

        if quiet is False:
            status_msg = '[{0}] athena annotate-pairs: Adding sequence homology ' + \
                         'features from reference fasta "{1}".'
            print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                    fasta))

        for identity in homology_cutoffs:
            for rev in True, False:
                pairs_df['newtrack_{}'.format(track_counter)] = \
                    add_homology(pairs_bt, fasta, binsize, identity, rev)
                track_counter += 1
    

    # Clean up long floats
    pairs_df = float_cleanup(pairs_df, maxfloat, start_idx=3)


    # Return bins as pbt.BedTool
    return pbt.BedTool.from_dataframe(pairs_df)

