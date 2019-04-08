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
import pybedtools


# Open connection to UCSC MySQL database for a specified reference
def ucsc_connect(build):
    conv = { MySQLdb.FIELD_TYPE.LONG: int }
    db = _mysql.connect(host='genome-mysql.cse.ucsc.edu', user='genome', 
                        passwd='', db=build, conv=conv)
    return db


# Check if a table exists
def table_exists(db, table):
    db.query('SELECT count(*) FROM information_schema.TABLES WHERE TABLE_NAME = "{0}"'.format(table))
    check = int(db.store_result().fetch_row(how=0, maxrows=0)[0][0])
    if check == 0:
        return False
    else:
        return True


# Convert chroms in a bedtool to hg18/hg19 nomenclature
def _check_hg_compliance(feature):
        if 'chr' not in feature.chrom:
            feature.chrom = 'chr' + feature.chrom
        return feature


# Convert chroms in a bedtool to GRC nomenclature
def _check_grc_compliance(feature):
        if 'chr' in feature.chrom:
            feature.chrom = str(feature.chrom).replace('chr', '')
        return feature


# Collapse bins to minimal ucsc query intervals
def collapse_query_regions(bins):
    query_ranges = bins.cut(range(3)).sort().merge()
    query_ranges = query_ranges.each(_check_hg_compliance).saveas()
    return query_ranges


# Add regional restrictions to a UCSC query SELECT statement
def constrain_query_regions(query, query_ranges):

    def _write_single_constraint(interval):
        c = '(`chrom` = "{0}" AND `chromEnd` >= {1} AND `chromStart` <= {2})'
        c = c.format(interval.chrom, interval.start, interval.stop)
        return c

    constraint = ''
    cidx = 0
    for interval in query_ranges:
        newcons = _write_single_constraint(interval)
        if cidx == 0:
            constraint = constraint + newcons
        else:
            constraint = ' '.join([constraint, 'OR', newcons])
        cidx += 1

    query = ' '.join([query, 'WHERE', constraint])

    return query


# Format UCSC query dependent on desired output format
def compose_query(table, db, oformat, query_ranges=None, map_column=None):

    # Query dependent on desired output behavior
    if oformat == 'bed':

        if map_column is None:
            query = 'SELECT `chrom`, `chromStart`, `chromEnd` from `{0}`'.format(table)

        else:
            query = 'SELECT `chrom`, `chromStart`, `chromEnd`, `{0}` from `{1}`'.format(map_column, table)
        
        # Constraint query region, if optioned
        if query_ranges is not None:
            query = constrain_query_regions(query, query_ranges)

    else:
        query = 'SELECT * from `{0}`'.format(table)

    return query


# Download a single table from UCSC database
def query_table(table, db, oformat='bed', query_ranges=None, map_column=None):

    # Query database
    query = compose_query(table, db, oformat, query_ranges, map_column)
    db.query(query)

    # Format output
    if oformat == 'bed':
        result = pybedtools.BedTool(db.store_result().fetch_row(how=0, maxrows=0))
        
    elif oformat == 'bigwig':
        result = dict(db.store_result().fetch_row(how=1, maxrows=0)[0])
        result = 'http://hgdownload.soe.ucsc.edu' + result['fileName'].decode('utf-8')
        
    else:
        result = db.store_result().fetch_row(how=1, maxrows=0)

    return result


# Master function for handling UCSC queries
def query_ucsc(bins, track, db, action, query_ranges, map_column=None):

    # Get raw data
    if 'map-' in action:
        if map_column is None:
            result = query_table(track, db, 'bigwig', query_ranges)
        else:
            result = query_table(track, db, 'bed', query_ranges, map_column)
            # Coerce to GRC nomenclature, if necessary
            if 'chr' not in bins[0]:
                result = result.each(_check_grc_compliance)

    elif action in 'count count-unique coverage'.split():
        result = query_table(track, db, 'bed', query_ranges)
        # Coerce to GRC nomenclature, if necessary
        if 'chr' not in bins[0]:
            result = result.each(_check_grc_compliance)

    else:
        from sys import exit
        exit('ERROR: Unknown action "{0}" in ucsc.query_ucsc'.format(action))

    return result
