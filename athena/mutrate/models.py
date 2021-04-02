#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Model architectures for mutation rate prediction
"""


import torch


class LogitModel(torch.nn.Module):
    """
    Defines a logistic regression model in torch.nn syntax
    """

    def __init__(self, n_features):
        super(LogitModel, self).__init__()
        self.linear = torch.nn.Linear(n_features, 1)

    def forward(self, x):
        y_pred = torch.sigmoid(self.linear(x))
        return y_pred


def initialize_torch_model(model_class, features, params={}):
    """
    Initialize PyTorch model, optimizer, and loss criteria based on command-line input
    """

    if model_class == 'logit':
        model = LogitModel(features.shape[1]).float()
        optimizer = torch.optim.SGD(model.parameters(), 
                                    lr=params.get('lr', 0.01),
                                    weight_decay=params.get('l2', 0))
        criterion = torch.nn.BCELoss()

    return model, optimizer, criterion


def train_torch_model(features, labels, model, optimizer, criterion, epochs=10e6, 
                      stop_early=False, earlyStopping={}, seed=2021):
    """
    Trains a PyTorch model from tensors of features and labels
    """

    loss_over_time = []
    test_loss_over_time = []

    # Check elements of earlyStopping dict
    train_eps = earlyStopping.get('train_eps', 10e-8)
    if 'features' in earlyStopping.keys() \
    and 'labels' in earlyStopping.keys():
        check_test_set = True
        test_features = earlyStopping['features']
        test_labels = earlyStopping['labels']
        monitor = earlyStopping.get('monitor', 5)
        test_eps = earlyStopping.get('train_eps', 10e-8)
        overfit_at_epoch = epochs
    else:
        check_test_set = False

    torch.manual_seed(seed)

    for epoch in range(int(epochs)):
        model.train()
        optimizer.zero_grad()

        # Forward pass
        preds = model(features)

        # Compute loss
        loss = criterion(preds, labels)
        loss_over_time.append(loss.float().item())

        # Backward pass
        loss.backward()
        optimizer.step()

        # Log test stats for early stopping, if optioned
        if stop_early and check_test_set and epoch % monitor == 0:
            model.eval()
            with torch.no_grad():
                test_loss = criterion(model(test_features), test_labels)
                test_loss_over_time.append(test_loss.float().item())
                if len(test_loss_over_time) > 1:
                    if test_loss_over_time[-2] < test_loss_over_time[-1]:
                        overfit_at_epoch = epoch
            model.train()

        # Check early stopping criteria
        if stop_early:
            if len(loss_over_time) > 1:
                # Stop if marginal training loss has converged below train_eps
                dLoss_train = loss_over_time[-2] - loss_over_time[-1]
                if abs(dLoss_train) < train_eps:
                    stop_reason = 'train_eps'
                    break
            if check_test_set:
                # Stop if marginal testing loss has converged below test_eps
                if len(test_loss_over_time) > 1:
                    dLoss_test = test_loss_over_time[-2] - test_loss_over_time[-1]
                    if abs(dLoss_test) < test_eps:
                        stop_reason = 'test_eps'
                        break
                    # Stop if testing loss started to increase in a previous 
                    # epoch and has continued to increase
                    if dLoss_test < 0 \
                    and overfit_at_epoch < epoch \
                    and epoch % monitor == 0:
                        stop_reason = 'overfit'
                        break

    if epoch + 1 == epochs:
        stop_reason = 'max_epochs'

    # Compile training info and 
    training_info = {'epochs_trained' : epoch + 1,
                     'loss_over_time' : loss_over_time,
                     'test_loss_over_time' : test_loss_over_time,
                     'stop_reason' : stop_reason}

    return model, training_info
