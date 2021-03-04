#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Annotate pairs of bins
"""


from athena.utils.misc import calc_binsize, determine_filetype
from athena.utils import nuc
from athena.utils.dfutils import float_cleanup
from gzip import GzipFile
import pybedtools as pbt
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


def add_homology(pairs_bt, fasta, binsize, homology_cutoffs=[1.0, 0.9, 0.7, 0.5]):
    """
    Compute homology features for all pairs from a reference fasta
    """

    # Process pairs in serial
    for pair in pairs_bt:

        # Only keep pairs representing distinct bins
        if pair.length > binsize:

            # Get nucleotide sequences corresponding to both bins in the pair
            bins_bt = _split_pairs(pair, binsize)
            seql, seqr = nuc.get_seqs_from_bt(bins_bt, fasta)

            # Get longest streches of sequences at varying levels of homology
            fwd_bp = [nuc.longest_subseq_min_identity(seql, seqr, i) for i in homology_cutoffs]
            rev_bp = [nuc.longest_subseq_min_identity(seql, seqr[::-1], i) for i in homology_cutoffs]
            import pdb; pdb.set_trace()

        # TODO: add case handling for overlapping bins in pair
        else:
            pass




    import pdb; pdb.set_trace()


def annotate_pairs(pairs, chroms, ranges, fasta, binsize, homology_cutoffs, maxfloat, quiet):
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


    # Annotate pairs based on nucleotide content, if optioned
    if fasta is not None:

        if quiet is False:
            status_msg = '[{0}] athena annotate-pairs: Adding sequence homology ' + \
                         'features from reference fasta "{1}".'
            print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                    fasta))

        add_homology(pairs_bt, fasta, binsize, homology_cutoffs)
    

    # Clean up long floats
    pairs_df = float_cleanup(pairs_df, maxfloat, start_idx=3)


    # Return bins as pbt.BedTool
    return pbt.BedTool.from_dataframe(pairs_df)

