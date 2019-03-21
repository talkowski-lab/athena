#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Athena: a toolkit for exploring structural variation mutation rates and dosage sensitivity
"""


import click
from .utils import commands as utils_commands
from .mutrate import commands as mutrate_commands
from .dosage import commands as dosage_commands


# Utility functions
@click.group(name='General utilities')
def utilscli():
    """
    General utilities
    """
    pass

utilscli.add_command(utils_commands.filtervcf)
utilscli.add_command(utils_commands.vcfstats)
utilscli.add_command(utils_commands.makebins)
utilscli.add_command(utils_commands.countsv)

# Mutation rate modeling functions
@click.group(name='Mutation rate modeling')
def mutratecli():
    """
    Mutation rate modeling
    """
    pass

mutratecli.add_command(mutrate_commands.annotatebins)
mutratecli.add_command(mutrate_commands.annodecomp)

# Dosage sensitivity modeling functions
@click.group(name='Dosage sensitivity modeling')
def dosagecli():
    """
    Dosage sensitivity modeling
    """
    pass

dosagecli.add_command(dosage_commands.mutquery)


# Top-level click command group for cli
cli = click.CommandCollection(sources=[utilscli, mutratecli, dosagecli])


# Main block
if __name__ == '__main__':
    cli()