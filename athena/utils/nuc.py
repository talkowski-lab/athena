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
from numpy import ceil, mean
from Bio import pairwise2
from Bio.Seq import Seq


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
            seqs.append(seq.rstrip())
    
    if return_headers:
        return seqs, headers
    else:
        return seqs


def _dice_sequence(seq, window_size, step_size):
    """
    Dices a string into sequential substrings of $window_size in $step_size steps
    """

    N = len(seq)

    if N <= window_size:
        return [seq]

    else:
        return [seq[k:k + window_size] for k in range(0, N - window_size, step_size)]


def pairwise_identity(seq1, seq2):
    """
    Compute strict (gapless) pairwise identity between two bins
    """

    N = min([len(seq1), len(seq2)])

    matches = pairwise2.align.globalms(Seq(seq1), Seq(seq2), 1, 0, -10e10, -10e10, 
                                       one_alignment_only=True, score_only=True)

    return matches / N


def sample_kmer_identity(seq1, seq2, k=100):
    """
    Compute the average strict (gapless) pairwise itentity between sequential k-mers of seq2 vs. seq1
    """

    kmers = _dice_sequence(seq2, k, k)

    identities = []
    for kmer in kmers:
        matches = pairwise2.align.localms(Seq(kmer), Seq(seq1), 1, 0, -10e10, -10e10, 
                                          one_alignment_only=True, score_only=True)
        identities.append(matches / k)

    return mean(identities), max(identities)


def longest_subseq_min_identity(seq1, seq2, min_identity=1.0, window_size=100):
    """
    Finds longest stretch of contiguous (gapless) sequences with at least $identity homology
    """

    ref = Seq(seq1)
    # step_size = int(ceil(window_size / 2))

    # # Step across seq2 to count total number of matches
    # seeds = _dice_sequence(seq2, window_size, step_size)
    # seed_matches = [pairwise2.align.localms(ref, Seq(query), 1, 0, -10e10, -10e10, score_only=True, one_alignment_only=True) for query in seeds]

    # test10 = [x for x in set(_dice_sequence(seq2, 10, 1)) if x in seq1]
    # test50 = [x for x in set(_dice_sequence(seq2, 50, 1)) if x in seq1]
    # test100 = [x for x in set(_dice_sequence(seq2, 100, 1)) if x in seq1]

    import pdb; pdb.set_trace()
    alignments = pairwise2.align.globalms(Seq(seq1), Seq(seq2), 1, 0, -10e10, -10e10, one_alignment_only=True, score_only=True)

    




