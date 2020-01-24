#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Read and sanitize a subset of columns from a .tsv (e.g., annotations from a BED)
"""

import warnings
import numpy as np
import pandas as pd


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


def load_feature_df(tsv, first_column=3, log_transform=None, sqrt_transform=None,
                    fill_missing=None, warn=False):
    """
    Wrapper function to load and sanitize features from a tsv input
    """

    # Disable data scaling warnings
    if not warn:
        warnings.filterwarnings(action='ignore', category=DataConversionWarning)

    # Read all columns from BED file except for 0:first_column
    ncols = len(pd.read_csv(tsv, sep='\t', nrows=1).columns)
    df = pd.read_csv(tsv, sep='\t', usecols=range(first_column, ncols))

    # Clean missing values
    if fill_missing is not None:
        df = _clean_missing_vals(df, fill_missing)

    # Transform any columns as optioned
    if log_transform is not None:
        for column in log_transform:
            df = _log_trans(df, column)

    if sqrt_transform is not None:
        for column in sqrt_transform:
            df = _sqrt_trans(df, column)

    return df

