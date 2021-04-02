#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Train a mutation rate model
"""


from athena.utils import dfutils
from athena.mutrate import models
import torch
import torch.utils.data as data_utils
import csv
from datetime import datetime
import random


def load_training_bed(bed_path, as_tensor=True):
    """
    Load a single training bed as a tuple of Tensors (features, labels)
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


def load_all_training_beds(pairs_tsv, as_tensor=True):
    """
    Wrapper to load all training BED files
    """

    data_dict = {}

    with open(pairs_tsv) as fin:
        for contig, bed_path in csv.reader(fin, delimiter='\t'):
            features, labels = load_training_bed(bed_path, as_tensor)
            data_dict[contig] = {'features' : features, 'labels' : labels}

    return data_dict


def make_training_tensors(data_dict, xchroms=[], seed=2021):
    """
    Pool torch.tensors from multiple chromosomes for model training
    """

    feats_list = [vals['features'] for chrom, vals in data_dict.items() if chrom not in xchroms]
    labs_list = [vals['labels'] for chrom, vals in data_dict.items() if chrom not in xchroms]

    train_features = torch.cat(feats_list, 0)
    train_labels = torch.cat(labs_list, 0)

    # Must shuffle features and labels to avoid nonuniformities during optimization
    torch.manual_seed(seed)
    new_order = torch.randperm(train_labels.shape[0])
    train_features = train_features[new_order, :]
    train_labels = train_labels[new_order, :]

    return train_features, train_labels


def eval_model(pred, labels, cutoff=0.5, noise=10e-10):
    """
    Evaluates performance of trained model predictions vs. real labels
    """

    pred_yes = (pred >= cutoff)
    label_yes = (labels >= cutoff)

    TP = (pred_yes & label_yes).sum().item()
    FP = (pred_yes & ~label_yes).sum().item()
    TN = (~pred_yes & ~label_yes).sum().item()
    FN = (~pred_yes & label_yes).sum().item()
    P = max([TP + FN, noise])
    N = max([TN + FP, noise])

    acc = (TP + TN) / (P + N)
    sens = TP / P
    spec = TN / N
    PPV = TP / max([(TP + FP), noise])
    NPV = TN / max([(TN + FN), noise])
    F1 = (2 * PPV * sens) / max([(PPV + sens), noise])
    
    stats = {'accuracy' : acc,
             'sensitivity' : sens,
             'specificity' : spec,
             'PPV' : PPV,
             'NPV' : NPV,
             'F1' : F1}

    return stats


def contig_cv(data_dict, test_contig, model_class, seed=2021):
    """
    Wrapper to conduct a single round of training & testing for a held-out chromosome
    """

    # Split data
    train_contigs = [k for k in data_dict.keys() if k != test_contig]
    train_features, train_labels = \
        make_training_tensors(data_dict, xchroms=test_contig, seed=seed)
    test_features, test_labels = \
        make_training_tensors(data_dict, xchroms=train_contigs, seed=seed)

    # Fit model
    model, optimizer, criterion = \
        models.initialize_torch_model(model_class, train_features)
    fit_model = models.train_torch_model(train_features, train_labels, model, 
                                         optimizer, criterion, seed=seed)

    # # Compute training & testing statistics
    model.eval()
    with torch.no_grad():
        train_stats = eval_model(fit_model(train_features), train_labels)
        test_stats = eval_model(fit_model(test_features), test_labels)
    model.train()

    return model, train_stats, test_stats


def mu_train(training_data, model_class, outfile, stats_outfile, cv_eval, 
             max_cv_k, quiet, seed=2021):
    """
    Master function to train a single mutation rate model
    """

    # Load and process all training BEDs
    if not quiet:
        status_msg = '[{0}] athena mu-predict: Loading training datasets.'
        print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S')))
    data_dict = load_all_training_beds(training_data)

    # Assign chromosomes to be held out for cross-validation
    cv_k = min([len(data_dict), max_cv_k])
    random.seed(seed)
    cv_test_contigs = sorted(random.sample(data_dict.keys(), cv_k))
    if not quiet:
        status_msg = '[{0}] athena mu-predict: Holding out data for the following ' + \
                     '{1} contigs as cross-validation test sets: {2}'
        print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'),
                                cv_k, ', '.join(cv_test_contigs)))

    # Evaluate training & testing performance with cross-validation
    for test_contig in cv_test_contigs:
        cv_res = {}
        fit_model, train_stats, test_stats = \
            contig_cv(data_dict, test_contig, model_class, seed=seed)
        cv_res[test_contig] = {'model' : fit_model,
                               'train_stats' : train_stats,
                               'test_stats' : test_stats}
        if not quiet:
            status_msg = '[{0}] athena mu-predict: Cross-validation results for {1}:'
            print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                    test_contig))
            print('  TRAINING:')
            for key, value in train_stats.items():
                print('    * {} = {}'.format(key, round(value, 6)))
            print('  TESTING:')
            for key, value in test_stats.items():
                print('    * {} = {}'.format(key, round(value, 6)))

    # After cross-validation, train model on all data for best predictive power
    all_features, all_labels = \
        make_training_tensors(data_dict, xchroms=[], seed=seed)
    model, optimizer, criterion = \
        models.initialize_torch_model(model_class, all_features)
    final_model = models.train_torch_model(all_features, all_labels, model, 
                                           optimizer, criterion)

    # Evaluate training performance of final model
    final_model.eval()
    final_stats = eval_model(final_model(all_features), all_labels)
    if not quiet:
        status_msg = '[{0}] athena mu-predict: Training performance for full model:'
        print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S')))
        for key, value in final_stats.items():
            print('    * {} = {}'.format(key, round(value, 6)))

    # Save trained model to outfile (note: model is intentionally exported in .eval() mode)
    torch.save(final_model, outfile)

    # Compile & save training stats, if optioned
