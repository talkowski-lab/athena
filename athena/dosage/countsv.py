#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019-Present Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Count overlap of SVs and generic BED or GTF files
"""


from athena.mutrate.countsv import _load_sv_from_bed
from athena.mutrate.muquery import _get_query_entity_name
from athena.utils.misc import determine_filetype, vcf2bed
import pybedtools as pbt
import pysam
import pandas as pd


def count_sv_generic(sv_in, query_in, outfile, group_by, ovr_frac, maxfloat, bgzip):
    """
    Collect SV counts for every feature in query_in 
    """

    # Load query file as pbt.BedTool
    qbt = pbt.BedTool(query_in)

    # Load SVs and convert to pbt.BedTool
    sv_format = determine_filetype(sv_in)
    if 'vcf' in sv_format:
        vcf = pysam.VariantFile(sv_in)
        sv = vcf2bed(vcf, breakpoints=False, add_ci_to_bkpts=False)
    elif 'bed' in sv_format:
        sv = _load_sv_from_bed(sv_in, breakpoints=False)

    # Iterate over all hits in intersection of query & SVs
    # For each hit, write SV ID to dict per query feature
    hit_ids = {}
    for hit in qbt.intersect(sv, loj=True, f=ovr_frac):
        qname  = _get_query_entity_name(hit, group_by)
        if qname not in hit_ids.keys():
            hit_ids[qname] = set()
        svid = hit[-1]
        if svid != '.':
            hit_ids[qname].add(svid)

    # Once hits have been processed, write out counts as .tsv
    hit_counts = {qname: len(svids) for qname, svids in hit_ids.items()}
    hit_df = pd.DataFrame.from_dict(hit_counts, orient='index').reset_index()
    hit_df.columns = '#query n_svs'.split()
    hit_df.to_csv(outfile, index=False, header=True, sep='\t')

