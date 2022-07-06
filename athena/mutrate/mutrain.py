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
from athena.mutrate import mututils
import torch
import torch.utils.data as data_utils
from datetime import datetime
import random
from numpy import median
import pandas as pd


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
    
    stats = {'total_positives' : TP + FN,
             'total_negatives' : TN + FP, 
             'accuracy' : acc,
             'sensitivity' : sens,
             'specificity' : spec,
             'PPV' : PPV,
             'NPV' : NPV,
             'F1' : F1}

    return stats


def contig_cv(data_dict, test_contig, model_class, hypers):
    """
    Wrapper to conduct a single round of training & testing for a held-out chromosome
    """

    seed = hypers['seed']

    # Split data
    train_contigs = [k for k in data_dict.keys() if k != test_contig]
    train_features, train_labels = \
        mututils.pool_tensors(data_dict, xchroms=test_contig, seed=seed)
    test_features, test_labels = \
        mututils.pool_tensors(data_dict, xchroms=train_contigs, seed=seed)

    # Fit model with early stopping
    model, optimizer, criterion = \
        models.initialize_torch_model(model_class, train_features, hypers)
    earlyStopping = {'features' : test_features, 'labels' : test_labels, 'monitor' : 5}
    fit_model, training_info = \
        models.train_torch_model(train_features, train_labels, model, optimizer, 
                                 criterion, stop_early=True, 
                                 earlyStopping=earlyStopping, seed=seed)

    # # Compute training & testing statistics
    model.eval()
    with torch.no_grad():
        train_stats = eval_model(fit_model(train_features), train_labels)
        test_stats = eval_model(fit_model(test_features), test_labels)
    model.train()

    return model, training_info, train_stats, test_stats


def interpret_stop_reason(stop_reason):
    """
    Convert 'stop_reason' code into plaintext explanation
    """

    if stop_reason == 'train_eps':
        stop_explain = 'reaching minimum change in training loss'
    elif stop_reason == 'test_eps':
        stop_explain = 'reaching minimum change in testing loss'
    elif stop_reason == 'overfit':
        stop_explain = 'testing loss beginning to increase'
    elif stop_reason == 'max_epochs':
        stop_explain = 'reaching maximum number of epochs'

    return stop_explain


def make_stats_df(cv_res, final_stats, avg_epochs):
    """
    Make a table containing all training & testing stats from cross-validation
    Also includes training stats for full model
    """

    cols = '#test_chrom stage epochs stop_reason total_positives total_negatives ' + \
           'accuracy sensitivity specificity ppv npv f1'
    cols = cols.split()
    stats_df = pd.DataFrame(columns=cols)

    # Add train & test stats for each round of CV
    for contig in cv_res.keys():
        for stage in 'train test'.split():
            vals = cv_res[contig][stage + '_stats']
            row_vals = [contig, stage, cv_res[contig]['epochs'], \
                        cv_res[contig]['stop_reason']] + list(vals.values())
            newrow = pd.Series(row_vals, index=cols)
            stats_df = stats_df.append(newrow, ignore_index=True)
    
    # Add training stats for full model
    row_vals = ['None', 'train', avg_epochs, 'max_epochs'] + \
               list(final_stats.values())
    newrow = pd.Series(row_vals, index=cols)
    stats_df = stats_df.append(newrow, ignore_index=True)
    
    return stats_df


def make_calibration_df(model, all_features, all_labels):
    """
    Apply trained model to training bins to generate calibration data
    """

    model.eval()
    with torch.no_grad():
        preds = model(all_features)
    model.train()
    
    cal_df = pd.DataFrame(torch.cat([preds, all_labels], 1).numpy(), 
                          columns='#predicted actual'.split())
    
    return cal_df


def mu_train(training_data, model_class, model_out, stats_out, cal_out, hypers, 
             maxfloat, quiet):
    """
    Master function to train a single mutation rate model
    """

    # Unpack certain hypers
    seed = hypers['seed']

    # Load and process all training BEDs
    if not quiet:
        status_msg = '[{0}] athena mu-train: Loading training datasets.'
        print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S')))
    data_dict = mututils.load_all_beds(training_data)

    # Perform cross-validation, if optioned
    if hypers['cv_eval']:

        # Assign chromosomes to be held out for cross-validation
        cv_k = min([len(data_dict), hypers['max_cv_k']])
        random.seed(seed)
        cv_test_contigs = sorted(random.sample(data_dict.keys(), cv_k))
        if not quiet:
            status_msg = '[{0}] athena mu-train: Holding out data for the following ' + \
                         '{1} contigs as cross-validation test sets: {2}'
            print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'),
                                    cv_k, ', '.join(cv_test_contigs)))

        # Evaluate training & testing performance with cross-validation
        cv_res = {}
        for test_contig in cv_test_contigs:
            fit_model, training_info, train_stats, test_stats = \
                contig_cv(data_dict, test_contig, model_class, hypers)
            stop_reason = training_info['stop_reason']
            cv_res[test_contig] = {'model' : fit_model,
                                   'train_stats' : train_stats,
                                   'test_stats' : test_stats,
                                   'epochs' : training_info['epochs_trained'],
                                   'stop_reason' : stop_reason}
            if not quiet:
                stop_explain = interpret_stop_reason(stop_reason)
                status_msg = '[{0}] athena mu-train: Cross-validation results for {1} ' + \
                             'after {2} epochs (stopped due to {3}):'
                print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                        test_contig, training_info['epochs_trained'], stop_explain))
                print('  TRAINING:')
                for key, value in train_stats.items():
                    print('    * {} = {}'.format(key, round(value, 6)))
                print('  TESTING:')
                for key, value in test_stats.items():
                    print('    * {} = {}'.format(key, round(value, 6)))

    # After cross-validation, train model on all data for best predictive power
    # Stop after median number of epochs from CV, above
    all_features, all_labels = \
        mututils.pool_tensors(data_dict, xchroms=[], seed=seed)
    model, optimizer, criterion = \
        models.initialize_torch_model(model_class, all_features, hypers)
    if hypers['cv_eval']:
        avg_epochs = int(median([vals['epochs'] for vals in cv_res.values()]))
    else:
        avg_epochs = hypers.get('max_epochs', 10e6)
    final_model, training_info = \
        models.train_torch_model(all_features, all_labels, model, optimizer, 
                                 criterion, stop_early=False, epochs=avg_epochs)

    # Evaluate training performance of final model
    final_model.eval()
    final_stats = eval_model(final_model(all_features), all_labels)
    if not quiet:
        status_msg = '[{0}] athena mu-train: Training performance for full model ' + \
                     'after {1} epochs:'
        print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), avg_epochs))
        for key, value in final_stats.items():
            print('    * {} = {}'.format(key, round(value, 6)))

    # Save trained model to model_out, if optioned
    # (note: model is intentionally exported in .eval() mode)
    if model_out is not None:
        torch.save(final_model, model_out)

    # Compile & save training stats, if optioned
    if stats_out  is not None:
        stats_df = make_stats_df(cv_res, final_stats, avg_epochs)
        stats_df = dfutils.float_cleanup(stats_df, maxfloat=maxfloat, start_idx=6)
        stats_df.to_csv(stats_out, sep='\t', index=False)

    # Compile & save data for calibration analysis, if optioned
    if cal_out is not None:
        cal_df = make_calibration_df(model, all_features, all_labels)
        cal_df = dfutils.float_cleanup(cal_df, maxfloat=maxfloat, start_idx=0)
        cal_df.to_csv(cal_out, sep='\t', index=False)
