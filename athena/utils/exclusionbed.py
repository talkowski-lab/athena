#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021- Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Utilities for handling exclusion BED files
"""

import pybedtools as pbt


def _buffer_exclusion_list(interval, excl_buffer=0):
    interval.start = max([0, interval.start - excl_buffer])
    interval.stop = interval.stop + excl_buffer
    return interval


def load_exclusion_bts(exclusion_list, excl_buffer=0):
    """
    Merge and buffer all exclusion lists into a single pbt.BedTool
    """

    n_xl = len(exclusion_list)

    if n_xl == 0:
        xl = pbt.BedTool('', from_string=True).saveas()
    else:
        xl = pbt.BedTool(exclusion_list[0])
        if n_xl > 1:
            xl = xl.cat(*[pbt.BedTool(x) for x in exclusion_list[1:]])

    return xl.each(_buffer_exclusion_list, excl_buffer=excl_buffer).merge().saveas()
