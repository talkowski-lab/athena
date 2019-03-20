#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""


import click


@click.command(name='query')
def mutquery():
    """
    Mutation rate lookup
    """
    click.echo('Query mutation rate matrix by interval (in dev.)')
