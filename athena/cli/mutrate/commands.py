#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""


import click


@click.command(name='annotate')
def annotatebins():
    """
    Annotate bins
    """
    click.echo('Annotate bins with one or more external tracks (in dev.)')

@click.command(name='decomp')
def annodecomp():
    """
    Decompose bin annotations
    """
    click.echo('Perform eigenvector decomposition on annotated bins (in dev.)')
