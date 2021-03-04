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


def longest_subseq_min_identity(seq1, seq2, min_identity):
    """
    Finds longest stretch of contiguous (gapless) sequences with at least $identity homology
    """

    import pdb; pdb.set_trace()
    alignments = pairwise2.align.localms(Seq(seq1), Seq(seq2), 
                                         1,                 # match score
                                         -100*min_identity, # mismatch penalty
                                         -10e10,            # gap open score
                                         -10e10,            # gap extend score
                                          one_alignment_only=True)

    




