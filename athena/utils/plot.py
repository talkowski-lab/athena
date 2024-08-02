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
        plot_title = title.replace('/', '_').replace(' ', '_')
        _simple_hist(vals, title)
        plt.savefig('.'.join([png_prefix, plot_title, 'png']), format='png')


def feature_importance(pca, feature_names, pdf_prefix, norm_variance=False,
                       abs_val=False, sort_features=False, sort_pcs=False,
                       pc_weights=None):
    """
    Plot matrix of raw feature importances through PCs
    """

    # Get matrix of raw feature weights in PCs
    weights = pca.components_

    # If optioned, normalize feature weights by variance explained by respective PCs
    if norm_variance:
        weights = weights * pca.explained_variance_[:, None]

    # If optioned, normalize feature weights by additional specified PC weights
    if pc_weights is not None:
        n_pc_weights = pc_weights.shape[0]
        if n_pc_weights < weights.shape[0]:
            # Truncate feature weights to specified number of PCs
            weights = weights[:n_pc_weights]
        weights = weights * pc_weights[:, None]

    n_pcs = weights.shape[0]

    # Prepare to plot additional panel with PC weights/variance explained
    if pc_weights is not None:
        bar_weights = pc_weights
        bar_title = "Weights on PCs"
    elif norm_variance:
        bar_weights = pca.explained_variance_
        bar_title = "PC % variance explained"
    else:
        bar_weights = None
        bar_title = None

    # If optioned, take absolute value for ease of interpretability
    if abs_val:
        weights = abs(weights)

    # Set heatmap colorscale midpoint at 0
    if abs_val:
        colors = "Purples"
        col_min = 0
        col_max = np.max(weights)
    else:
        colors = "bwr"
        col_min = -np.max(abs(weights))
        col_max = np.max(abs(weights))

    # If optioned, reorder rows and columns according to weights
    feature_labels = feature_names
    pc_labels = np.arange(n_pcs) + 1
    if sort_features:
        col_idx = weights.max(axis=0).argsort()[::-1]
        weights = weights[:, col_idx]
        feature_labels = [feature_labels[i] for i in col_idx]
    if sort_pcs:
        row_idx = weights.max(axis=1).argsort()[::-1]
        weights = weights[row_idx, :]
        pc_labels = [pc_labels[i] for i in row_idx]

    # Set plot dimensions
    plot_cols = 2 if bar_weights is not None else 1
    width_ratio = [3, 1] if bar_weights is not None else [1]
    fig, axs = plt.subplots(
        1, plot_cols, sharey=True, layout="constrained",
        gridspec_kw={'width_ratios': width_ratio},
        squeeze=False
    )
    axs = axs.reshape(-1)
    fig.set_figheight(10)

    # Plot weight matrix
    axs[0].imshow(weights, cmap=colors, vmin=col_min, vmax=col_max)

    # Replace x-labels with feature names
    axs[0].set_xticks(np.arange(len(feature_labels)), labels=feature_labels)
    axs[0].tick_params(axis="x", labelrotation=90)
    axs[0].set_yticks(np.arange(len(pc_labels)), labels=pc_labels)

    # Add axes & title
    axs[0].set_xlabel("Raw features")
    axs[0].set_ylabel("PCs")
    title = "Raw feature weights in PCs"
    if norm_variance:
        title = title + "\n* PC variance explained"
    if pc_weights is not None:
        title = title + "\n* weights on PCs"
    axs[0].set_title(title)

    # Add additional panel for PC weights/explained variance
    if bar_weights is not None:
        if sort_pcs:
            bar_weights = bar_weights[row_idx]
        # Reverse bar weight order as they are plotted from bottom to top
        axs[1].barh(
            np.arange(len(pc_labels))[::-1], bar_weights[::-1],
            align="center"
        )
        # Add axes & title
        axs[1].set_title(bar_title)
        axs[1].tick_params(axis="y", length=0)
        # Resize to match image size
        asp = abs(
            (np.diff(axs[1].get_xlim())[0] / np.diff(axs[1].get_ylim())[0])
            / (np.diff(axs[0].get_xlim())[0] / np.diff(axs[0].get_ylim())[0])
            / (width_ratio[1] / width_ratio[0])
        )
        axs[1].set_aspect(asp)

    plt.savefig(
        ".".join([pdf_prefix, "feature_importance", "pdf"]),
        format="pdf",
        bbox_inches="tight",
    )
