#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Apply a pre-trained mutation rate model to predict mu for new bins
"""


import pandas as pd
from athena.utils import dfutils
from athena.utils.misc import determine_filetype, bgzip as bgz
from athena.mutrate.mututils import load_bed
import torch
from numpy import log10
from os import path


def mu_predict(pairs, model_pkl, outfile, raw_mu, keep_features, maxfloat, bgzip):
    """
    Apply a trained mutation rate model to new bin-pairs
    """

    # Load pairs and split coordinates from features
    coords = pd.read_csv(pairs, sep='\t', usecols=range(3))
    feats, labels = load_bed(pairs)
    if keep_features:
        feats_df = dfutils.load_feature_df(pairs)

    # Load model from .pkl and switch to evaluation mode
    model = torch.load(model_pkl)
    model.eval()

    # Predict mutation rates for all bins
    with torch.no_grad():
        preds = model(feats).numpy()
        if not raw_mu:
            preds = log10(preds)
        preds_df = pd.DataFrame(preds, columns=['mu'])

    # Format output dataframe
    out_df = pd.concat([coords, preds_df], axis=1)
    if keep_features:
        out_df = pd.concat([out_df, feats_df], axis=1)
    out_df = dfutils.float_cleanup(out_df, maxfloat=maxfloat, start_idx=3)

    # Save pairs with predicted mutation rates
    if 'compressed' in determine_filetype(outfile):
        outfile = path.splitext(outfile)[0]
    out_df.to_csv(outfile, sep='\t', index=False)

    # Bgzip bins, if optioned
    if bgzip:
        bgz(outfile)
