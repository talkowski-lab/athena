#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Perform Eigendecomposition on 1D bin annotations
"""


import pandas as pd
import numpy as np
from os import path
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA, FastICA
from scipy.stats import spearmanr
import pybedtools
from athena.utils.misc import bgzip as bgz
import athena.utils.dfutils as dfutils


# Compute PCA & feature stats
def get_feature_stats(df_annos, feature_names, pca, pcs, stats_outfile, 
                      eigen_prefix):
    # Compute variance explained
    expvar = list(pca.explained_variance_ratio_)

    # Prep feature correlation dictionary
    components = pca.n_components_
    eigen_names = ['_'.join([eigen_prefix, str(i+1)]) for i in range(components)]
    corstats = {}
    for name in eigen_names:
        if name not in corstats.keys():
            corstats[name] = []

    # Determine correlations of each raw annotation vs each PC
    for comp in range(components):
        ename = eigen_names[comp]
        for feat in range(df_annos.shape[1]):
            fname = feature_names[feat]
            if np.ptp(df_annos, axis=0)[feat] > 0:
                spear = tuple(spearmanr(pcs[:, comp], df_annos[:, feat]))
                corstats[ename].append((fname, spear[0], spear[1]))
            else:
                corstats[ename].append((fname, np.nan, np.nan))

    # Format results and write to text file
    fout = open(stats_outfile, 'w')

    for comp in range(components):

        header_line = '## Eigenfeature {0} ({1}% variance explained)\n'
        fout.write(header_line.format(comp + 1, round(100 * expvar[comp], 2)))

        ename = eigen_names[comp]

        for fname, rho, pval in sorted(corstats[ename], key=lambda x: (x[-1], -abs(x[1]))):
            if pval <= 0.05 / len(feature_names):
                newline = '    * {0}: Rho = {1}; P = {2:.2e}\n'
                fout.write(newline.format(fname, round(rho, 2), pval))

        fout.write('\n\n')


# Clean up long floats in a data frame
def float_cleanup(df, maxfloat, start_idx):
    df.iloc[:, start_idx:] = df.iloc[:, start_idx:].replace('.', np.nan)
    df.iloc[:, start_idx:] = df.iloc[:, start_idx:].apply(pd.to_numeric)
    df.iloc[:, start_idx:] = df.iloc[:, start_idx:].round(maxfloat)
    return df


# Master function for Eigendecomposition of bin annotations
def decompose_bins(bins, outfile, ica=False, components=10, minvar=None, 
                   log_transform=None, sqrt_transform=None, exp_transform=None, 
                   square_transform=None, boxcox_transform=None, fill_missing=0, 
                   first_column=3, maxfloat=5, max_pcs=100, pca_stats=None, 
                   eigen_prefix='eigenfeature', bgzip=False):

    # Read bins, then sanitize and transform annotations
    df_bins = pd.read_csv(bins, sep='\t', usecols=range(first_column))
    df_annos = dfutils.load_feature_df(bins, first_column, log_transform,
                                       sqrt_transform, exp_transform, 
                                       square_transform,  boxcox_transform,
                                       fill_missing)
    feature_names = df_annos.columns.tolist()

    # Scale all columns
    df_annos = StandardScaler().fit_transform(df_annos)

    # Decompose annotations
    pcs_to_calc = min([df_annos.shape[1], max_pcs])
    if ica:
        ica = FastICA(n_components=pcs_to_calc)
        ics = ica.fit_transform(df_annos)
        exit('DEV NOTE: ICA is not currently supported. Exiting.')
    else:
        pca = PCA(n_components=pcs_to_calc)
        pcs = pca.fit_transform(df_annos)
        if minvar is None:
            components = pcs_to_calc
        else:
            components = len([i for i in np.cumsum(pca.explained_variance_ratio_) \
                              if i < minvar])

    eigen_names = ['_'.join([eigen_prefix, str(i+1)]) for i in range(components)]
    df_pcs = pd.DataFrame(pcs[:, :components], columns = eigen_names)

    # Write output bins with PCs
    if '.gz' in outfile:
        outfile = path.splitext(outfile)[0]
    out_df = float_cleanup(pd.concat([df_bins, df_pcs], axis=1), maxfloat, first_column)
    out_df.to_csv(outfile, sep='\t', index=False)
    if bgzip:
        bgz(outfile)

    # Perform extra assessments of PCA & feature fits, if optioned
    if pca_stats is not None:
        get_feature_stats(df_annos, feature_names, pca, pcs, pca_stats, eigen_prefix)

