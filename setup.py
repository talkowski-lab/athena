#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Setup script for Athena
"""

from setuptools import setup, find_packages

setup(
    name='athena',
    version='0.1dev',
    license='MIT',
    description='Athena: a toolkit for exploring structural variation ' + 
                'mutation rates and dosage sensitivity',
    long_description=open('README.md').read(),
    author='Ryan L. Collins',
    author_email='rlcollins@g.harvard.edu',
    packages=['athena'],
    include_package_data=True,
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        athena=athena.cli.master:cli
    ''',
)