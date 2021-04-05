#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Various utilities for mutation rate modeling
"""

from athena.utils import dfutils
import torch
import csv


def load_bed(bed_path, as_tensor=True):
    """
    Load a single BED of bins with SV counts and features as a tuple of Tensors (features, labels)
    """

    df = dfutils.load_feature_df(bed_path)

    if 'sv' not in df.columns:
        err = 'ERROR: required column "sv" not found in input BED {}'
        from sys import exit
        exit(err.format(bed_path))

    if as_tensor:
        features = torch.tensor(df.drop('sv', axis=1).values).float()
        labels = torch.tensor(df[['sv']].values).float()
        return features, labels


def load_all_beds(pairs_tsv, as_tensor=True):
    """
    Wraps load_bed() to import all BED files for a list of chromosomes
    """

    data_dict = {}

    with open(pairs_tsv) as fin:
        for contig, bed_path in csv.reader(fin, delimiter='\t'):
            features, labels = load_bed(bed_path, as_tensor)
            data_dict[contig] = {'features' : features, 'labels' : labels}

    return data_dict


def pool_tensors(data_dict, xchroms=[], seed=2021):
    """
    Pool torch.tensors from multiple chromosomes for model training or application
    """

    feats_list = [vals['features'] for chrom, vals in data_dict.items() if chrom not in xchroms]
    labs_list = [vals['labels'] for chrom, vals in data_dict.items() if chrom not in xchroms]

    all_features = torch.cat(feats_list, 0)
    all_labels = torch.cat(labs_list, 0)

    # Must shuffle features and labels to avoid nonuniformities during optimization
    torch.manual_seed(seed)
    new_order = torch.randperm(all_labels.shape[0])
    all_features = all_features[new_order, :]
    all_labels = all_labels[new_order, :]

    return all_features, all_labels

