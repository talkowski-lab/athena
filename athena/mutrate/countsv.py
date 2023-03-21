#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019-Present Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Count overlap of SVs and 1D bins or 2D bin-pairs
"""


import pybedtools as pbt
import pysam
from os import path
import gzip
from athena.utils.misc import bgzip as bgz
from athena.utils.misc import determine_filetype, calc_binsize, vcf2bed, add_names_to_bed
from athena.utils import dfutils
import pandas as pd
import numpy as np
from scipy.stats import norm
from sys import stdout


def _load_sv_from_bed(sv_in, breakpoints=False):
    """
    Loads SVs from an input BED file and optionally splits into breakpoints
    """

    intervals = ''
    bkpt_fmt = '{}\t{}\t{}\t{}\t{}\t{}\n'

    for interval in pbt.BedTool(sv_in).cut(range(3)):
        chrom = interval.chrom
        start = interval.start
        end = interval.end
        vid = '{}_{}_{}'.format(chrom, start, end)

        if breakpoints:
            left_interval = bkpt_fmt.format(chrom, start, start + 1, vid, 'POS', start)
            right_interval = bkpt_fmt.format(chrom, end, end + 1, vid, 'END', end)
            new_interval = left_interval + right_interval

        else:
            new_interval = '{}\t{}\t{}\t{}\n'.format(chrom, start, end, vid)

        intervals += new_interval

    return pbt.BedTool(intervals, from_string=True)


def _split_pairs(pairs_bt, binsize, add_name=False, add_side=False):
    """
    Splits a pbt.BedTool of bin-pairs into its constituent bins
    """

    pairs_str = ''

    for rec in pairs_bt:

        chrom, left, right = rec.fields[0:3]
    
        lbin = '\t'.join([chrom, left, str(int(left) + binsize)])
        rbin = '\t'.join([chrom, str(int(right) - binsize), right])

        if add_name:
            lbin = '\t'.join([lbin, rec.name])
            rbin = '\t'.join([rbin, rec.name])

        if add_side:
            lbin += '\tLEFT'
            rbin += '\tRIGHT'

        lbin += '\n'
        rbin += '\n'

        pairs_str += lbin + rbin

    return pbt.BedTool(pairs_str, from_string=True)


def calc_prob_bkpt(b_start, b_end, sv_mid, sv_start_sd, sv_end_sd):
    """
    Apportions breakpoint probability density across the interval b_start to b_end
    """

    b_start_centered = b_start - sv_mid
    b_end_centered = b_end - sv_mid

    # Note: left and right tails of distribution must be treated separately because
    # SV breakpoint uncertainty is not always symmetric (and can be encoded in VCF as such)

    if sv_start_sd > 0:
        left_z_min = min([-100, b_start_centered / sv_start_sd])
        left_z_max = min([0, b_end_centered / sv_end_sd])
        left_tail_p = norm.cdf(left_z_max) - norm.cdf(left_z_min)
    else:
        left_tail_p = 0.5

    if sv_end_sd > 0:
        right_z_min = max([0, b_start_centered / sv_start_sd])
        right_z_max = min([100, max([0, b_end_centered / sv_end_sd])])
        right_tail_p = norm.cdf(right_z_max) - norm.cdf(right_z_min)
    else:
        right_tail_p = 0.5

    return left_tail_p + right_tail_p


def _resolve_2d_bkpt_probs(p_dict):
    """
    Helper function to handle the combination of breakpoint probabilities for a
    single comparison of 2D bin-pairs
    """

    left = p_dict.get('LEFT', {})
    right = p_dict.get('RIGHT', {})
    lr_svs = set(left.keys()).intersection(set(right.keys()))

    # Compute the combined probability for each SV
    p_per_sv = np.array([])
    for sv_id in lr_svs:
        posend = p_dict['LEFT'][sv_id]['POS'] * p_dict['RIGHT'][sv_id]['END']
        endpos = p_dict['LEFT'][sv_id]['END'] * p_dict['RIGHT'][sv_id]['POS']
        sv_p = 1 - np.prod(1 - np.array([posend, endpos]))
        p_per_sv = np.append(p_per_sv, sv_p)

    # Compute the overall probability of _at least one_ SV existing between the bin-pair
    return 1 - np.prod(1 - p_per_sv)


def parse_breakpoint_hits(hits, paired, probs):
    """
    Summarize breakpoint counts or probabilities per bin from a precomputed intersection
    """

    bkpt_res_tmp = {}

    # Process each hit one at a time
    for hit in hits:

        # Define variables from hit
        b_chrom, b_start, b_end, b_id = hit.fields[:4]
        if paired:
            b_side = hit.fields[4]
            sv_chrom, sv_start, sv_end, sv_id = hit.fields[5:9]
            sv_side, sv_mid, sv_start_sd, sv_end_sd, z_extend = hit.fields[9:]
        else:
            sv_chrom, sv_start, sv_end, sv_id = hit.fields[4:8]
            sv_side, sv_mid, sv_start_sd, sv_end_sd, z_extend = hit.fields[8:]

        # If reporting breakpoint probabilities, need to store marginal probability per breakpoint per bin
        if probs:
            calc_prob_bkpt_inputs = [float(x) for x in [b_start, b_end, sv_mid, sv_start_sd, sv_end_sd]]
            prob_bkpt = calc_prob_bkpt(*calc_prob_bkpt_inputs)
            if paired:
                if b_id not in bkpt_res_tmp.keys():
                    bkpt_res_tmp[b_id] = {}
                if b_side not in bkpt_res_tmp[b_id].keys():
                    bkpt_res_tmp[b_id][b_side] = {}
                if sv_id not in bkpt_res_tmp[b_id][b_side]:
                    bkpt_res_tmp[b_id][b_side][sv_id] = {'POS' : 0, 'END' : 0}
                bkpt_res_tmp[b_id][b_side][sv_id][sv_side] = prob_bkpt
            else:
                bkpt_res_tmp[b_id] = np.append(bkpt_res_tmp.get(b_id, np.array([])), prob_bkpt)

        # If reporting breakpoint counts, add sv_id to set() of breakpoints per bin
        else:
            if paired:
                if b_id not in bkpt_res_tmp.keys():
                    bkpt_res_tmp[b_id] = {}
                if b_side not in bkpt_res_tmp[b_id].keys():
                    bkpt_res_tmp[b_id][b_side] = {}
                if sv_side not in bkpt_res_tmp[b_id][b_side].keys():
                    bkpt_res_tmp[b_id][b_side][sv_side] = set()
                bkpt_res_tmp[b_id][b_side][sv_side].add(sv_id)    

            else:
                if b_id not in bkpt_res_tmp.keys():
                    bkpt_res_tmp[b_id] = set()
                bkpt_res_tmp[b_id].add(sv_id)

    # Consolidate results
    if probs:
        if paired:
            # For 2D bins, need to first compute the sum of LEFT/RIGHT and POS/END 
            # probability products per SV per bin, then compute the complement of 
            # the product of the complement of _those_ probabilities
            # This is done in a helper function for clarity:
            bkpt_res = {b_id : _resolve_2d_bkpt_probs(p_dict) for b_id, p_dict in bkpt_res_tmp.items()}

        else:
            # For 1D bins, if multiple breakpoints overlap the same bin,
            # compute the complement of the product of their probability complements
            # i.e., the probability that strictly zero of the breakpoints truly
            # lie within the bin
            bkpt_res = {b_id : round(1 - np.prod(1 - prob_arr), 6) for b_id, prob_arr in bkpt_res_tmp.items()}

    else:
        if paired:
            bkpt_res = {}
            for b_id, sv_hits in bkpt_res_tmp.items():
                left_pos = sv_hits.get('LEFT', {}).get('POS', set())
                left_end = sv_hits.get('LEFT', {}).get('END', set())
                right_pos = sv_hits.get('RIGHT', {}).get('POS', set())
                right_end = sv_hits.get('RIGHT', {}).get('END', set())
                lr_posend = left_pos.intersection(right_end)
                lr_endpos = left_end.intersection(right_pos)
                bkpt_res[b_id] = len(lr_posend.union(lr_endpos))

        else:
            bkpt_res = {b_id : len(sv_ids) for b_id, sv_ids in bkpt_res_tmp.items()}

    return bkpt_res


def count_sv_in_bins(sv_in, bins_in, outfile, paired, binsize, breakpoints, 
                     probs, sv_ci, maxfloat, bgzip):
    """
    Main function to annotate bins_in with count (or probability) of SVs
    """

    # Load bins, split bin coordinates from annotations, and retain header
    if 'compressed' in determine_filetype(bins_in):
        bins_header = gzip.open(bins_in, 'r').readline().decode('utf-8').rstrip().split('\t')
    else:
        bins_header = open(bins_in, 'r').readline().rstrip().split('\t')
    bins_bt = pbt.BedTool(bins_in).cut(range(3)).saveas()
    bins_df = bins_bt.to_dataframe()
    feats_df = dfutils.load_feature_df(bins_in)
    if binsize is None:
        binsize = calc_binsize(bins_in)

    # Parse input SV file depending on format
    # If breakpoints == False, will return simple four-column BED with variant ID in fourth column
    # If breakpoints == True, will return two rows per record where each record
    # is one breakpoint with columns 4 = variant ID, 5 = POS or END, 6 = original
    # POS or END coordinate, 7 = std dev of left side of breakpoint, 8 = std dev of
    # right side of breakpoint, and 9 = number of std deviations extended left & right (i.e., z_extend)
    sv_format = determine_filetype(sv_in)
    if 'vcf' in sv_format:
        vcf = pysam.VariantFile(sv_in)
        sv = vcf2bed(vcf, breakpoints=breakpoints, add_ci_to_bkpts=probs, ci=sv_ci)
    elif 'bed' in sv_format:
        sv = _load_sv_from_bed(sv_in, breakpoints=breakpoints)

    # Perform intersection with bins depending on input parameters
    if breakpoints:
        bins_bt = add_names_to_bed(bins_bt)
        bin_ids = [b.name for b in bins_bt]

        # Split pairs if necessary
        if paired:
            bins_bt = _split_pairs(bins_bt, binsize=binsize, add_name=True, add_side=True)

        # Intersect breakpoints with bins
        hits = bins_bt.intersect(sv, wa=True, wb=True)
        bkpt_res = parse_breakpoint_hits(hits, paired, probs)
        sv_column = pd.Series([bkpt_res.get(b_id, 0) for b_id in bin_ids])

    # --comparison "overlap" (i.e., breakpoints == False) is the same for both 1D and 2D bins
    else:
        if probs:
            sv_column = pd.Series([min([1, int(x[-1])]) for x in bins_bt.intersect(sv, c=True)])
        else:
            sv_column = pd.Series([int(x[-1]) for x in bins_bt.intersect(sv, c=True)])

    # Paste bin coordinates, SV counts, and original features into single dataframe
    out_df = dfutils.float_cleanup(pd.concat([bins_df, sv_column, feats_df], axis=1), 
                                   maxfloat, 3)
    out_df.columns = bins_header[:3] + ['sv'] + bins_header[3:]

    # Save bins with SV counts
    if outfile in 'stdout /dev/stdout -'.split():
        outfile = stdout
        outfile_is_stdout = True
    else:
        outfile_is_stdout = False
        if 'compressed' in determine_filetype(outfile):
            outfile = path.splitext(outfile)[0]
    out_df.to_csv(outfile, sep='\t', header=True, index=False)

    # Bgzip bins, if optioned
    if bgzip and not outfile_is_stdout:
        bgz(outfile)

