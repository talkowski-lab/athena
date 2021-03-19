#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021- Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Utility functions for extracting and manipulating nucleotide sequences
"""


from athena.utils.misc import determine_filetype
from gzip import GzipFile
import pybedtools as pbt
import itertools
from rapidfuzz import fuzz


def get_seqs_from_bt(bt, fasta, return_headers=False):
    """
    Extract a list of nucleotide sequences corresponding to all records in a pbt.BedTool
    """

    if 'compressed' in determine_filetype(fasta):
        fasta = GzipFile(fasta)
    fseqs = bt.sequence(fasta).seqfn

    seqs = []
    headers = []
    with open(fseqs) as fin:
        for seqheader, seq in itertools.zip_longest(*[fin]*2):
            headers.append(seqheader.rstrip().replace('>', ''))
            seqs.append(seq.rstrip().upper())
    
    if return_headers:
        return seqs, headers
    else:
        return seqs


def longest_subseq_min_identity(seq1, seq2, min_identity=1.0):
    """
    Finds longest stretch of contiguous (gapless) sequences with at least $identity homology
    """

    seq1_len = len(seq1)
    longest = 0
    start = 0
    end = start + longest + 1
    fuzz_cutoff = 100 * min_identity

    while end <= seq1_len + 1:
        
        query = seq2[start:end]
        
        if min_identity == 1:
            match = query in seq1
        else:
            match = fuzz.partial_ratio(query, seq1) >= fuzz_cutoff

        if match:
            longest = len(query)
            end += 1
        else:
            start += 1

        end += 1
    
    return longest
    
