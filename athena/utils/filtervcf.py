#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Filter an input VCF
"""


import pysam
from sys import stdin, stdout, exit
from os import path
import pybedtools
import numpy as np
from athena.utils.misc import bgzip as bgz
from athena.utils.misc import hwe_chisq


def filter_vcf(vcf, out, chroms, xchroms, svtypes, exclusion_list, 
               minAF, maxAF, minAC, maxAC, minAN, filters, 
               minQUAL, maxQUAL, HWE, af_field, keep_infos, bgzip):

    # Open connection to input VCF
    if vcf in '- stdin'.split():
        invcf = pysam.VariantFile(stdin) 
    else:
        invcf = pysam.VariantFile(vcf)
    header = invcf.header

    #Clean undesired INFO fields from header
    if keep_infos != 'ALL':
        if keep_infos is None:
            keep_infos = []
        else:
            keep_infos = keep_infos.split(',')
        for key in 'END CHR2 SVTYPE SVLEN'.split():
            keep_infos.append(key)
        for key in header.info.keys():
            if key not in keep_infos:
                header.info.remove_header(key)

    # Open connection to output VCF
    if out in '- stdout'.split():
        outvcf = pysam.VariantFile(stdout, 'w', header=header)
    else:
        if '.gz' in out:
            out = path.splitext(out)[0]
        outvcf = pysam.VariantFile(out, 'w', header=header)

    # Parse filtering options
    if chroms is not None:
        chroms = chroms.split(',')
    else:
        chroms = header.contigs.keys()
    if xchroms is not None:
        xchroms = xchroms.split(',')
        chroms = [c for c in chroms if c not in xchroms]
    if svtypes is not None:
        if 'SVTYPE' not in header.info.keys():
            sys.exit('SVTYPE filtering was specified, but input VCF ' +
                     'does not have SVTYPE entry in INFO.')
        else: 
            svtypes = svtypes.split(',')
    if filters is not None:
        filters = filters.split(',')
    if exclusion_list is not None:
        bl = pybedtools.BedTool(exclusion_list)
        # bl = pd.read_csv(exclusion_list, sep='\t', names='chrom start end'.split())

    # Raise warning if AF or AC are missing from VCF
    for key in [af_field, 'AC']:
        if key not in header.info.keys():
            import warnings
            warning_message = '{0} not found in VCF INFO, so {0}-based filtering ' + \
                              'will be ignored'
            warning_message = warning_message.format(key)
            warnings.warn(warning_message, RuntimeWarning)
    
    # Raise exception if HWE enabled but any necessary fields missing
    if HWE is not None:
        for key in 'N_HOMREF N_HET N_HOMALT'.split():
            if key not in header.info.keys():
                error_message = 'Hardy-Weinberg filtering not possible due to ' + \
                                'missing {0} in VCF INFO'
                sys.exit(error_message.format(key))

    # Iterate over vcf & filter records
    for record in invcf.fetch():
        # Filter by chromosome
        if chroms is not None \
        and record.chrom not in chroms:
            continue

        # Filter by svtype
        if svtypes is not None \
        and record.info['SVTYPE'] not in svtypes:
            continue

        # Exclude records where end < start
        if record.stop < record.start:
            continue

        # Filter by AF/AC
        if af_field in record.info.keys():
            if minAF is not None:
                if np.nansum(record.info[af_field]) < minAF:
                    continue
            if maxAF is not None:
                if np.nansum(record.info[af_field]) > maxAF:
                    continue        
        if 'AC' in record.info.keys():
            if minAC is not None:
                if np.nansum(record.info['AC']) < minAC:
                    continue
            if maxAC is not None:
                if np.nansum(record.info['AC']) > maxAC:
                    continue

        # Filter by AN
        if 'AN' in record.info.keys():
            if minAN is not None:
                if record.info['AN'] < minAN:
                    continue

        # Filter by VCF FILTER
        if filters is not None:
            if len([f for f in record.filter if f not in filters]) > 0:
                continue

        # Filter by QUAL score
        if minQUAL is not None \
        and record.qual is not None:
            if record.qual < minQUAL:
                continue
        if maxQUAL is not None \
        and record.qual is not None:
            if record.qual > maxQUAL:
                continue

        # Filter by Hardy-Weinberg equilibrium
        if HWE is not None and len(record.alts) < 3:
            if np.nansum(record.info[af_field]) < 1:
                if hwe_chisq(record) < HWE:
                    continue

        # Clean record
        if keep_infos != 'ALL':
            for key in record.info.keys():
                if key not in keep_infos:
                    record.info.pop(key)

        # Write filter-passing records to output VCF
        outvcf.write(record)

    outvcf.close()

    # Filter remaining records against exclusion_list
    if exclusion_list is not None:
        prebl_vcf = pybedtools.BedTool(out)
        prebl_vcf.intersect(bl, header=True, v=True).saveas(out)

    # Bgzip output VCF, if optioned
    if bgzip:
        bgz(out)

