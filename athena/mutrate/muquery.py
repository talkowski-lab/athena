#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Query a mutation rate matrix
"""


from os import path
from athena.utils.misc import determine_filetype, check_contig_naming_scheme
from athena.utils.dfutils import float_cleanup
import re
import pybedtools as pbt
import pysam
import pandas as pd
from numpy import log10, nansum, nanmax
from sys import stdout


def _get_query_entity_name(feature, query_group_by):
    """
    Determine query entity name from pbt.Interval
    """

    if query_group_by is not None:
        query_name = feature.attrs[query_group_by]
    else:
        if 'gene_name' in feature.attrs.keys():
            query_name = feature.attrs['gene_name']
        elif len(feature.fields) > 3:
            query_name = feature.name
        else:
            query_name = '_'.join([feature.chrom, str(feature.start), str(feature.end)])

    return query_name


def mu_query(pairs, query, outfile, query_group_by, ovr_frac, raw_mu_in, 
             raw_mu_out, epsilon, maxfloat, gzip):
    """
    Query a mutation rate matrix
    """

    # Infer contig naming scheme of mutrates from head of pairs
    mu_contig_naming = check_contig_naming_scheme(pairs)

    # Preprocess query
    if path.exists(query):
        query_ftype = determine_filetype(query)
        if query_ftype in 'bed gtf compressed-bed compressed-gtf'.split():
            qbt = pbt.BedTool(query)
        else:
            err = 'Determined file type "{}" for query file {}; this file type ' + \
                  'is currently not supported in athena mu-query.'
            from sys import exit
            exit(err.format(query_ftype, query))

    else:
        query = re.sub('[\t\:-]', '|', query).split('|')
        qbt = pbt.BedTool('\t'.join(query) + '\n', from_string=True)

    # Iterate over qbt and extract mutation rates for each
    # Builds results as dict of pd.DataFrame of slices of mutation rate matrix
    # keyed on entity name (or coordinates if no names are found)
    query_results_dfs = {}
    qres_columns ='chrom start end mu'.split()
    with pysam.TabixFile(pairs) as mutrates:
        for qint in qbt:

            # Check for contig naming consistency between query & mutrates
            if mu_contig_naming == 'no_chr' and 'chr' in qint.chrom:
                qint.chrom = re.sub('^chr', '', str(qint.chrom))
            elif mu_contig_naming == 'has_chr' and 'chr' not in qint.chrom:
                qint.chrom = 'chr' + str(qint.chrom)

            # Determine query entity name
            query_name = _get_query_entity_name(qint, query_group_by)
            if query_name not in query_results_dfs.keys():
                query_results_dfs[query_name] = pd.DataFrame(columns=qres_columns)

            # Query mutation rate matrix
            qstr = '{}:{}-{}'.format(qint.chrom, qint.start, qint.end)
            qhits = [i for i in mutrates.fetch(qstr)]

            # Enforce minimum overlap fraction, if optioned
            if ovr_frac is not None:
                qres_bt = pbt.BedTool('\n'.join(qhits), from_string=True)
                qint_bt = pbt.BedTool(re.sub('[\:-]', '\t', qstr) + '\n', from_string=True)
                qres_hits = qres_bt.intersect(qint_bt, F=ovr_frac, wa=True, u=True)
                qres = qres_hits.to_dataframe(names=qres_columns)

            else:
                qres = [i.rstrip().split('\t') for i in qhits]
                qres = pd.DataFrame(qres, columns=qres_columns)    

            # Add passing mutation rate cells to dict of results dataframes
            query_results_dfs[query_name] = \
                pd.concat([query_results_dfs[query_name], qres], ignore_index=True)

    # Deduplicate results for each query entity and sum mutation rates
    query_results = {}
    for qname, mu_df in query_results_dfs.items():
        keepers = ~mu_df['chrom start end'.split()].duplicated()
        mu_df = mu_df.loc[keepers, :].copy()
        if not raw_mu_in:
            mu_df.loc[:, 'mu'] = 10 ** mu_df['mu'].astype(float)
        mu_sum = nanmax([epsilon, nansum(mu_df['mu'])])
        if not raw_mu_out:
            mu_sum = log10(mu_sum)
        query_results[qname] = mu_sum
    query_results = pd.DataFrame.from_dict(query_results, orient='index').reset_index()
    query_results.columns = '#query mu'.split()
    query_results = float_cleanup(query_results, maxfloat, 1)
    
    # Write query to output file
    if outfile in 'stdout /dev/stdout -'.split():
        outfile = stdout
    query_results.to_csv(outfile, header=True, index=False, sep='\t')

