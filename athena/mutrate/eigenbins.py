#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019-Present Ryan L. Collins <rlcollins@g.harvard.edu>
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
import pickle
import pybedtools
from athena.utils.misc import bgzip as bgz
from athena.utils.misc import determine_filetype
import athena.utils.dfutils as dfutils


def _load_precomp_model(precomp_model):
    """
    Load a precomputed eigen-bins model from a .pickle to be applied to new data
    """

    with open(precomp_model, 'rb') as pkl_in:
        df_fills, trans_dict, scaler, pca, components, whitener, eigenval_limits \
            = pickle.load(pkl_in)

    return df_fills, trans_dict, scaler, pca, components, whitener, eigenval_limits


def _save_model_params(df_fills, trans_dict, scaler, pca, components, whitener, 
                       eigenval_limits, parameters_outfile):
    """
    Save a model's parameters as a .pickle for application to other data
    """

    with open(parameters_outfile, 'wb') as pkl_out:
        model_to_save = [df_fills, trans_dict, scaler, pca, components, 
                         whitener, eigenval_limits]
        pickle.dump(model_to_save, pkl_out)


def get_feature_stats(df_annos, feature_names, pca, pcs, stats_outfile, 
                      eigen_prefix, components=None):
    """
    Compute PCA & feature stats
    """

    # Compute variance explained
    expvar = list(pca.explained_variance_ratio_)

    # Prep feature correlation dictionary
    if components is None:
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


def decompose_bins(bins, bins_outfile=None, parameters_outfile=None, precomp_model=None, 
                   components=10, minvar=None, trans_dict=None, whiten=False, 
                   cap_predictions=False, fill_missing=0, first_column=3, 
                   maxfloat=5, max_pcs=100, pca_stats=None, 
                   eigen_prefix='eigenfeature', bgzip=False):
    """
    Master function for Eigendecomposition of bin annotations
    """

    # Set certain defaults prior to loading precomputed model
    whitener = None

    # Load precomputed model, if optioned
    if precomp_model is not None:
        df_fills, trans_dict, scaler, pca, components, whitener, eigenval_limits = \
            _load_precomp_model(precomp_model)
        fill_missing = df_fills

    # Expand feature transformation dictionary
    log_transform = trans_dict.get('log', [])
    sqrt_transform = trans_dict.get('sqrt', [])
    exp_transform = trans_dict.get('exp', [])
    square_transform = trans_dict.get('square', [])
    boxcox_transform = trans_dict.get('boxcox', [])

    # Read bins, then sanitize and transform annotations
    df_bins = pd.read_csv(bins, sep='\t', usecols=range(first_column))
    df_annos, df_fills = \
        dfutils.load_feature_df(bins, first_column, log_transform, sqrt_transform, 
                                exp_transform, square_transform,  boxcox_transform, 
                                fill_missing, return_fills=True)
    feature_names = df_annos.columns.tolist()
    
    # Scale all columns
    if precomp_model is None:
        scaler = StandardScaler().fit(df_annos)
    df_annos = scaler.transform(df_annos)

    # Learn covariance matrix & determine number of components to keep
    if precomp_model is None:
        pcs_to_calc = min([df_annos.shape[1], max_pcs])
        pca = PCA(n_components=pcs_to_calc).fit(df_annos)
        if minvar is None:
            components = pcs_to_calc
        else:
            components = len([i for i in np.cumsum(pca.explained_variance_ratio_) \
                              if i < minvar])

    # Decompose annotations
    pcs = pca.transform(df_annos)
    eigen_names = ['_'.join([eigen_prefix, str(i+1)]) for i in range(components)]
    df_pcs = pd.DataFrame(pcs[:, :components], columns=eigen_names)

    # "Whiten" eigenfeatures, if optioned
    if whiten:
        if precomp_model is None:
            whitener = StandardScaler().fit(df_pcs)
    if whitener is not None:
        df_pcs = pd.DataFrame(whitener.transform(df_pcs), columns=eigen_names)

    # Cap extreme eigenvalues, if optioned
    if precomp_model is None:
        smallest_eigenvals = pd.Series([-np.inf] * df_pcs.shape[1], index=df_pcs.columns)
        largest_eigenvals = pd.Series([np.inf] * df_pcs.shape[1], index=df_pcs.columns)
        if cap_predictions:
            smallest_eigenvals = df_pcs.min()
            largest_eigenvals = df_pcs.max()
        eigenval_limits = {'lower' : smallest_eigenvals,
                           'upper' : largest_eigenvals}
    dfutils.cap_values(df_pcs, eigenval_limits['lower'], direction='lower')
    dfutils.cap_values(df_pcs, eigenval_limits['upper'], direction='upper')

    # Write output bins with PCs
    if bins_outfile is not None:
        if 'compressed' in determine_filetype(bins_outfile):
            bins_outfile = path.splitext(bins_outfile)[0]
        out_df = dfutils.float_cleanup(pd.concat([df_bins, df_pcs], axis=1), 
                                       maxfloat, first_column)
        out_df.to_csv(bins_outfile, sep='\t', index=False)
        if bgzip:
            bgz(bins_outfile)

    # Save model for future use, if optioned
    if parameters_outfile is not None:
        _save_model_params(df_fills, trans_dict, scaler, pca, components, 
                           whitener, eigenval_limits, parameters_outfile)

    # Perform extra assessments of PCA & feature fits, if optioned
    if pca_stats is not None:
        get_feature_stats(df_annos, feature_names, pca, pcs, pca_stats, 
                          eigen_prefix, components)

