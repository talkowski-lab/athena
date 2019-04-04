#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Functions to interact with UCSC Genome Browser database
"""


import MySQLdb
import _mysql


# Open connection to UCSC MySQL database for a specified reference
def ucsc_connect(build):
    conv = { MySQLdb.FIELD_TYPE.LONG: int }
    db = _mysql.connect(host='genome-mysql.cse.ucsc.edu', user='genome', 
                        passwd='', db=build, conv=conv)
    return db


# Download a single table from UCSC database
def query_ucsc_table(table, db, format='bed'):
    
    # Behavior: convert table to simple BedTool
    if format == 'bed':
        query = 'SELECT chrom, chromStart, chromEnd from {0}'.format(table)
        db.query(query)
        result = pybedtools.BedTool(db.store_result().fetch_row(how=0, maxrows=0))

    return result




