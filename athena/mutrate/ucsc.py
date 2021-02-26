#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Functions to interact with UCSC Genome Browser database
"""


import MySQLdb
from MySQLdb import _mysql
import pybedtools
import re


default_columns = ['chrom', 'chromStart', 'chromEnd']


def ucsc_connect(build):
    """
    Open connection to UCSC MySQL database for a specified reference
    """
    conv = { MySQLdb.FIELD_TYPE.LONG: int }
    db = _mysql.connect(host='genome-mysql.cse.ucsc.edu', user='genome', 
                        passwd='', db=build, conv=conv)
    return db


def table_exists(db, table):
    """
    Check if a table exists
    """
    db.query('SELECT count(*) FROM information_schema.TABLES WHERE TABLE_NAME = "{0}"'.format(table))
    check = int(db.store_result().fetch_row(how=0, maxrows=0)[0][0])
    if check == 0:
        return False
    else:
        return True


def _check_hg_compliance(feature):
    """
    Convert chroms in a bedtool to hg18/hg19 nomenclature
    """

    if 'chr' not in feature.chrom:
        feature.chrom = 'chr' + feature.chrom
    return feature


def _check_grc_compliance(feature):
    """
    Convert chroms in a bedtool to GRC nomenclature
    """

    if 'chr' in feature.chrom:
        feature.chrom = str(feature.chrom).replace('chr', '')
    return feature


def collapse_query_regions(bins):
    """
    Collapse bins to minimal ucsc query intervals
    """

    query_ranges = bins.cut(range(3)).sort().merge()
    query_ranges = query_ranges.each(_check_hg_compliance).saveas()
    return query_ranges


def constrain_query_regions(query, columns, query_ranges):
    """
    Add regional restrictions to a UCSC query SELECT statement
    """

    def _write_single_constraint(interval):
        c = '(`{0}` = "{1}" AND `{2}` >= {3} AND `{4}` <= {5})'
        c = c.format(columns[0], interval.chrom,
                     columns[1], interval.start,
                     columns[2], interval.stop)
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

    if ' WHERE ' in query:
        query = '{0} AND ( {1} )'.format(query, constraint)
    else:
        query = ' '.join([query, 'WHERE', constraint])

    return query


def parse_table_arg(track):
    """
    Parse table & column input
    """

    if ':' in track:
        table = track.split(':')[0]
        track_opts = track.split(':')[1].split(',')
        
        columns = [i for i in track_opts if not any(s in i for s in '! = < >'.split())]
        if len(columns) < 3:
            columns = default_columns + columns

        conditions = [i for i in track_opts if any(s in i for s in '! = < >'.split())]
        if len(conditions) == 0:
            conditions = None

    else:
        table = track
        columns = default_columns
        conditions = None

    return table, columns, conditions


def format_conditions(conditions):
    """
    Format conditions into SQL-compliant syntax
    """

    def _form_cond(cond):
        terms = re.sub('[!=<>]', ' ', cond).split()
        col = '`{0}`'.format(terms[0])
        query = terms[1]
        if not query.isdigit():
            query = '"{0}"'.format(query)
        comp = ''.join([c for c in cond if c in set('!=><')])
        form_cond = ' '.join([col, comp, query])
        return form_cond

    form_conds = [_form_cond(cond) for cond in conditions]

    return form_conds


def compose_query(table, db, oformat, query_ranges=None, 
                  columns=default_columns, conditions=None):
    """
    Format UCSC query dependent on desired output format
    """

    # Query dependent on desired output behavior
    if oformat == 'bed':

        query = 'SELECT ' + ', '.join(['`{0}`'.format(c) for c in columns])

        query = query + ' from `{0}`'.format(table)
        
        # Add other conditions, as optioned
        if conditions is not None:
            form_conds = ' AND '.join(format_conditions(conditions))
            query = '{0} WHERE {1}'.format(query, form_conds)

        # Constraint query region, if optioned
        if query_ranges is not None:
            query = constrain_query_regions(query, columns, query_ranges)

    else:
        query = 'SELECT * from `{0}`'.format(table)

    return query


def query_table(table, db, oformat='bed', query_ranges=None, 
                columns=default_columns, conditions=None):
    """
    Download a single table from UCSC database
    """

    # Query database
    query = compose_query(table, db, oformat, query_ranges, columns, conditions)
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


def query_ucsc(bins, track, columns, conditions, db, action, query_ranges):
    """
    Master function for handling UCSC queries
    """

    # Get raw data
    if 'map-' in action:
        if columns is default_columns:
            result = query_table(track, db, 'bigwig', query_ranges)
        else:
            result = query_table(track, db, 'bed', query_ranges, 
                                 columns, conditions)
            # Coerce to GRC nomenclature, if necessary
            if 'chr' not in bins[0].chrom:
                result = result.each(_check_grc_compliance)

    elif action in 'count count-unique coverage'.split():
        result = query_table(track, db, 'bed', query_ranges,
                             columns, conditions)
        # Coerce to GRC nomenclature, if necessary
        if 'chr' not in bins[0].chrom:
            result = result.each(_check_grc_compliance).saveas()

    else:
        from sys import exit
        exit('ERROR: Unknown action "{0}" in ucsc.query_ucsc'.format(action))

    return result


def subquery_ucsc(bins, track, columns, conditions, db, action, query_ranges,
                  contig):
    """
    UCSC query restricted to a single contig
    """

    sub_ranges = query_ranges.filter(lambda f: f.chrom == contig).saveas()
    result = query_ucsc(bins, track, columns, conditions, db, action, sub_ranges)
    if isinstance(result, pybedtools.BedTool):
        result = result.saveas()

    return result

