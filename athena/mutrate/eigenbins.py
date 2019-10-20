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
import warnings
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.exceptions import DataConversionWarning
from scipy.stats import spearmanr
import pybedtools
from athena.utils.misc import bgzip as bgz


# Handle missing annotation values
def _clean_missing_vals(df, fill_missing):
    df = df.replace('.', np.nan)
    df = df.apply(pd.to_numeric, errors='coerce')
    if str(fill_missing).isdigit():
        df = df.fillna(value=fill_missing)
    elif fill_missing == 'mean':
        df = df.fillna(df.mean())
    elif fill_missing == 'median':
        df = df.fillna(df.median())
    else:
        from sys import exit
        err = 'ERROR: "{0}" not a recognized specification for --fill-missing'
        exit(err.format(fill_missing))

    return df


# Perform log-transformation of a single feature
def _log_trans(df, column):
    
    if column not in list(df.columns):
        from sys import exit
        err = 'ERROR: "{0}" is not a recognized feature name in input bins; ' + \
              'log transformation failed'
        exit(err.format(column))

    df[[column]] = df[[column]].apply(pd.to_numeric, errors='coerce')

    maxval = df[[column]].max().tolist()[0]

    df[[column]] = np.log10(df[[column]] + maxval/100)

    return df


# Perform square root-transformation of a single feature
def _sqrt_trans(df, column):
    
    if column not in list(df.columns):
        from sys import exit
        err = 'ERROR: "{0}" is not a recognized feature name in input bins; ' + \
              'log transformation failed'
        exit(err.format(column))

    df[[column]] = df[[column]].apply(pd.to_numeric, errors='coerce')

    df[[column]] = np.sqrt(df[[column]])

    return df


# Compute PCA & feature stats
def get_feature_stats(df_annos, feature_names, pca, pcs, stats_outfile):
    # Compute variance explained
    expvar = list(pca.explained_variance_ratio_)

    # Prep feature correlation dictionary
    components = pca.n_components_
    eigen_names = ['eigenfeature_{0}'.format(i+1) for i in range(components)]
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
def decompose_bins(bins, outfile, components=10, log_transform=None, 
                   sqrt_transform=None, fill_missing=0, first_column=3, 
                   maxfloat=5, pca_stats=None, bgzip=False):

    # Disable data scaling warnings
    warnings.filterwarnings(action='ignore', category=DataConversionWarning)

    # Read bins & subset to columns with annotations
    df_all = pd.read_csv(bins, sep='\t')
    df_bins = df_all.iloc[:, :first_column]
    df_annos = df_all.iloc[:, first_column:]
    feature_names = df_annos.columns.tolist()

    # Clean missing values
    df_annos = _clean_missing_vals(df_annos, fill_missing)

    # Transform any columns as optioned
    if log_transform is not None:
        for column in log_transform:
            df_annos = _log_trans(df_annos, column)

    if sqrt_transform is not None:
        for column in sqrt_transform:
            df_annos = _sqrt_trans(df_annos, column)

    # Scale all columns
    df_annos = StandardScaler().fit_transform(df_annos)

    # Perform PCA
    pca = PCA(n_components=components)
    pcs = pca.fit_transform(df_annos)
    eigen_names = ['eigenfeature_{0}'.format(i+1) for i in range(components)]
    df_pcs = pd.DataFrame(pcs, columns = eigen_names)

    # Write output bins with PCs
    if '.gz' in outfile:
        outfile = path.splitext(outfile)[0]
    out_df = float_cleanup(pd.concat([df_bins, df_pcs], axis=1), maxfloat, 3)
    out_df.to_csv(outfile, sep='\t', index=False)
    if bgzip:
        bgz(outfile)

    # Perform extra assessments of PCA & feature fits, if optioned
    if pca_stats is not None:
        get_feature_stats(df_annos, feature_names, pca, pcs, pca_stats)

