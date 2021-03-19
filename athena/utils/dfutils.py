#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019- Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Read and sanitize a subset of columns from a .tsv (e.g., annotations from a BED)
"""


import csv
import numpy as np
import pandas as pd
from scipy.stats import boxcox
import warnings
from sklearn.exceptions import DataConversionWarning
import pandas as pd
from os import path
from athena.utils import bgzip as bgz


def _load_transformations(trans_tsv):
    """
    Read all data transformations from a tsv file
    """

    all_trans = {}

    with open(trans_tsv) as fin:
        tsv = csv.reader(fin, delimiter='\t')
        for trans, feature in tsv:
            trans = trans.replace('--', '').split('-')[0]
            if trans in all_trans.keys():
                all_trans[trans].append(feature)
            else:
                all_trans[trans] = [feature]

    return all_trans


def _clean_missing_vals(df, fill_missing, return_fills=False):
    """
    Handle missing annotation values
    """
    df = df.replace('.', np.nan)
    df = df.apply(pd.to_numeric, errors='coerce')
    df_fills = None
    if fill_missing is not None:
        if isinstance(fill_missing, dict):
            df = df.fillna(fill_missing)
            df_fills = fill_missing
        elif str(fill_missing).isdigit():
            df = df.fillna(value=fill_missing)
            df_fills = {k : fill_missing for k in df.columns}
        elif fill_missing == 'mean':
            df_fills = df.mean().to_dict()
            df = df.fillna(df_fills)
        elif fill_missing == 'median':
            df_fills = df.median().to_dict()
            df = df.fillna(df_fills)
        else:
            from sys import exit
            err = 'ERROR: "{0}" not a recognized specification for --fill-missing'
            exit(err.format(fill_missing))

    if return_fills:
        return df, df_fills
    else:
        return df


def _log_trans(df, column, strict=False):
    """
    Perform log-transformation of a single feature
    """
    
    if column not in list(df.columns):
        err = 'ERROR: "{0}" is not a recognized feature name in input bins; ' + \
              'log transformation failed'
        if strict:
            from sys import exit
            exit(err.format(column))
        else:
            print(err.format(column))
            return df

    else:
        df[[column]] = df[[column]].apply(pd.to_numeric, errors='coerce')
        maxval = np.nanmax(df[[column]])
        df[[column]] = np.log10(df[[column]] + maxval/1000)
        return df


def _sqrt_trans(df, column, strict=False):
    """
    Perform square root-transformation of a single feature
    """
    
    if column not in list(df.columns):
        err = 'ERROR: "{0}" is not a recognized feature name in input bins; ' + \
              'square root transformation failed'
        if strict:
            from sys import exit
            exit(err.format(column))
        else:
            print(err.format(column))
            return df

    else:
        df[[column]] = df[[column]].apply(pd.to_numeric, errors='coerce')
        df[[column]] = np.sqrt(df[[column]])
        return df


def _exp_trans(df, column, strict=False):
    """
    Perform exponential-transformation of a single feature
    """
    
    if column not in list(df.columns):
        err = 'ERROR: "{0}" is not a recognized feature name in input bins; ' + \
              'exponential transformation failed'
        if strict:
            from sys import exit
            exit(err.format(column))
        else:
            print(err.format(column))
            return df

    else:
        df[[column]] = df[[column]].apply(pd.to_numeric, errors='coerce')
        df[[column]] = np.exp(df[[column]])
        return df


def _square_trans(df, column, strict=False):
    """
    Perform square-transformation of a single feature
    """
    
    if column not in list(df.columns):
        err = 'ERROR: "{0}" is not a recognized feature name in input bins; ' + \
              'square transformation failed'
        if strict:
            from sys import exit
            exit(err.format(column))
        else:
            print(err.format(column))
            return df

    else:
        df[[column]] = df[[column]].apply(pd.to_numeric, errors='coerce')
        df[[column]] = df[[column]] ** 2
        return df


def _boxcox_trans(df, column, strict=False):
    """
    Perform Box-Cox power transformation of a single feature
    """
    
    if column not in list(df.columns):
        err = 'ERROR: "{0}" is not a recognized feature name in input bins; ' + \
              'Box-Cox transformation failed'
        if strict:
            from sys import exit
            exit(err.format(column))
        else:
            print(err.format(column))
            return df

    else:
        df[[column]] = df[[column]].apply(pd.to_numeric, errors='coerce')
        maxval = np.nanmax(df[[column]])
        minval = np.nanmin(df[[column]])
        if minval >= 0:
            df[column] = boxcox(df[column] + maxval/1000)[0]
        else:
            err = 'ERROR: negative values of "{}" are not allowed in Box-Cox ' + \
                      'transformation'
            if strict:
                from sys import exit
                exit(err.format(column))
            else:
                print(err.format(column))

        return df


def float_cleanup(df, maxfloat, start_idx):
    """
    Clean up long floats in a subset of columns of bins
    """

    df.iloc[:, start_idx:] = \
        df.iloc[:, start_idx:].replace('.', np.nan).apply(pd.to_numeric).round(maxfloat)

    return df


def load_feature_df(tsv, first_column=3, log_transform=None, sqrt_transform=None,
                    exp_transform=None, square_transform=None, boxcox_transform=None, 
                    fill_missing=None, return_fills=False, warn=False, strict=False):
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
    df, df_fills = _clean_missing_vals(df, fill_missing, return_fills=True)

    # Transform any columns as optioned
    if log_transform is not None:
        for column in log_transform:
            df = _log_trans(df, column, strict)

    if sqrt_transform is not None:
        for column in sqrt_transform:
            df = _sqrt_trans(df, column, strict)

    if exp_transform is not None:
        for column in exp_transform:
            df = _exp_trans(df, column, strict)

    if square_transform is not None:
        for column in square_transform:
            df = _square_trans(df, column, strict)

    if boxcox_transform is not None:
        for column in boxcox_transform:
            df = _boxcox_trans(df, column, strict)

    if return_fills:
        return df, df_fills
    else:
        return df


def transform_df(bed_in, bed_out, first_column=3, log_transform=None, 
                sqrt_transform=None, exp_transform=None, square_transform=None, 
                boxcox_transform=None, bgzip=False, fill_missing=None, 
                warn=False, strict=False):
    """
    Wrapper function to load, transform, and write an annotated binset
    """

    df_bins = pd.read_csv(bed_in, sep='\t', usecols=range(first_column))
    df_annos = load_feature_df(bed_in, first_column, log_transform, sqrt_transform,
                               exp_transform, square_transform, boxcox_transform,
                               fill_missing, warn, strict)
    out_df = pd.concat([df_bins, df_annos], axis=1)

    if '.gz' in bed_out:
        bed_out = path.splitext(bed_out)[0]
    out_df.to_csv(bed_out, sep='\t', na_rep='NA', index=False)

    if bgzip:
        bgz(bed_out)

