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
        optimizer = torch.optim.SGD(model.parameters(), lr=params.get('lr', 0.01))
        criterion = torch.nn.BCELoss()

    return model, optimizer, criterion


def train_torch_model(features, labels, model, optimizer, criterion, epochs=1000, seed=2021):
    """
    Trains a PyTorch model from a TensorDataset of (features, labels)
    """

    torch.manual_seed(seed)

    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()

        # Forward pass
        preds = model(features)

        # Compute loss
        loss = criterion(preds, labels)

        # Backward pass
        loss.backward()
        optimizer.step()

    return model
