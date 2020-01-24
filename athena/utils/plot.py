#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Plotting functions embedded in Athena
"""


import athena.utils.dfutils as dfutils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def feature_hists(bed, png_prefix, skip_cols=3, log_transform=None, 
                   sqrt_transform=None, exp_transform=None, 
                   square_transform=None,  boxcox_transform=None,
                   fill_missing=None):
    """
    Plot simple histograms for all columns in a BED file
    """

    # Load & sanitize features
    df = dfutils.load_feature_df(bed, skip_cols, log_transform, sqrt_transform,
                                 exp_transform, square_transform, 
                                 boxcox_transform, fill_missing)
    
    def _simple_hist(vals, title):
        """
        Plot a single simple histogram of values
        """

        # Subset to values within the middle 99.9% of data
        n_nan = len(vals[np.isnan(vals)])
        vals_num = vals[~np.isnan(vals)]
        vlims = np.quantile(vals_num, q=[0.0005, 0.9995])
        n_outlier = len(vals_num[(vals_num < vlims[0]) | (vals_num > vlims[1])])
        vals_plot = vals_num[(vals_num >= vlims[0]) & (vals_num <= vlims[1])]

        # Plot & format histogram
        fig, ax = plt.subplots()
        n, bins, patches = plt.hist(vals_plot, 25)
        plt.subplots_adjust(top=0.8)
        
        # Add axes & title
        ax.set_xlabel(title)
        ax.set_ylabel('Bins')
        fulltitle = '\n'.join([title,
                               '{:,} total bins'.format(len(vals)),
                               '{:,} bins with missing values (not shown)'.format(n_nan),
                               '{:,} outlier bins (not shown)'.format(n_outlier)])
        ax.set_title(fulltitle)

    # Plot one histogram per column
    for i in range(len(df.columns)):
        title = df.columns[i]
        vals = df[title]
        _simple_hist(vals, title)
        plt.savefig('.'.join([png_prefix, title, 'png']), format='png')

