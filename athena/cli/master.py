#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Athena: a toolkit for exploring structural variation mutation rates and dosage sensitivity
"""


import click


# Top-level click command group for cli
@click.group()
def cli():
    pass

@click.command(name='filter-vcf')
def filtervcf():
    click.echo('Filter input vcf (in development)')

@click.command(name='make-bins')
def makebins():
    click.echo('Create bins (in development)')

cli.add_command(filtervcf)
cli.add_command(makebins)


# Main block
if __name__ == '__main__':
    cli()