#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Functions for interfacing with Google Cloud Storage objects within Athena
"""


import os
import subprocess
from datetime import datetime
from athena.utils.misc import determine_extension
from sys import stdout


def _find_remote_index(url, ftype):
    """
    Locate tabix or samtools bam/cram index for a remotely hosted gs:// URL
    """

    # Build expected index filename
    if ftype in 'compressed-vcf compressed-bed'.split():
        ipath = url + '.tbi'
    elif ftype in 'bam':
        ipath = url + '.bai'
    elif ftype in 'cram':
        ipath = url + crai
    else:
        err = 'INDEX LOCALIZATION FAILURE: unknown index format expected for ' + \
              'input file {}'
        exit(err.format(url))

    # Try to find index locally first (no need to download)
    local_ipath = os.path.basename(ipath)
    if os.path.exists(local_ipath):
        return local_ipath

    # If local copy of index can't be found, try to query gs:// path
    else:
        gs_ls_pipe = subprocess.Popen(['gsutil', 'ls', ipath], stdout=subprocess.PIPE)
        gspath = gs_ls_pipe.communicate()[0].decode('utf-8').rstrip()
        # If index is remote but not local, need to download the index locally
        if gs_ls_pipe.returncode == 0:
            status_msg = '[{0}] athena slice-remote: Found index for input ' + \
                         'file "{1}" hosted remotely but no local copy found. ' + \
                         'Downloading with gsutil now.'
            print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                                    url))
            gs_cp_ecode = os.system('gsutil cp {} ./'.format(ipath))
            if gs_cp_ecode != 0:
                err = 'INDEX LOCALIZATION FAILURE: gsutil cp failed for remotely ' + \
                      'hosted index {}'
                exit(err.format(ipath))
            return local_ipath

        else:
            err = 'INDEX LOCALIZATION FAILURE: unable to locate expected index ' + \
                  '{} for input file {}'
            exit(err.format(ipath, url))


def _parse_remote_inputs_tsv(urls_tsv, local_suffix='local_slice'):
    """
    Parse a tsv of remote files to localize
    """

    urls_dict = {}

    with open(urls_tsv) as f_in:
        for k, line in enumerate(f_in):

            url, mdata = line.rstrip().split('\t', 1)
            basename_orig = os.path.basename(url)
            ftype, fext = determine_extension(basename_orig, return_extension=True)
            sliced_ext = '{}.{}'.format(local_suffix, fext)
            basename_sliced = basename_orig[:-len(fext)] + sliced_ext
            index_path = _find_remote_index(url, ftype)
            
            # Confirm 
            if ftype not in 'compressed-bed compressed-vcf bam cram'.split():
                err = 'INPUT ERROR: format not recognized as tabix or samtools ' + \
                      'compliant for input file {}'
                exit(err.format(url))

            urls_dict[k] = {'url' : url,
                            'ftype' : ftype,
                            'metadata' : mdata,
                            'basename_orig' : basename_orig,
                            'basename_sliced' : basename_sliced,
                            'index_path' : index_path}
    
    return urls_dict


def _check_cram_fasta_input(urls_dict, ref_fasta):
    """
    Ensure reference FASTA file is provided if 
    """

    ftypes = [vals['ftype'] for vals in urls_dict.values()]
    if 'cram' in ftypes:
        if not os.isfile(ref_fasta):
            err = 'INPUT ERROR: input .tsv contains one or more CRAM files ' + \
                  'but --ref-fasta not specified'
            exit(err)


def _authenticate_remote_streaming():
    """
    Ensure $GCS_OAUTH_TOKEN is set (required for remote streaming with htslib)
    """

    if 'GCS_OAUTH_TOKEN' not in os.environ:
        err = 'GCS AUTHENTICATION ERROR: environment variable GCS_OAUTH_TOKEN ' + \
              'unassigned but is required for remote streaming.\nPlease set ' + \
              'GCS_OAUTH_TOKEN and try again.'
        exit(err)
    else:
        pass


def _localize_remote_slice(url_info, regions_bed, ref_fasta=None):
    """
    Localize a slice of a remote file using samtools or tabix
    """

    ftype = url_info['ftype']

    if ftype in 'cram bam'.split():
        slice_method = 'samtools'
    else:
        slice_method = 'tabix'

    if ftype == 'cram':
        loc_cmd = 'samtools view -C -M -T ' + ref_fasta + ' -L {} {} > {}'
        idx_cmd = 'samtools index -c {}'
    elif ftype == 'bam':
        loc_cmd = 'samtools view -b -M -L {} {} > {}'
        idx_cmd = 'samtools index -b {}'
    elif slice_method =='tabix':
        loc_cmd = 'tabix -h -T {} {} | bgzip -c > {}'
        idx_cmd = 'tabix -f {}'
    else:
        err = 'LOCALIZATION ERROR: unable to determine localization method for {}'
        exit(err.format(url_info['url']))
    loc_cmd = loc_cmd.format(regions_bed, url_info['url'], url_info['basename_sliced'])
    idx_cmd = idx_cmd.format(url_info['basename_sliced'])

    status_msg = '[{0}] athena slice-remote: Slicing {1} with {2}:\n{3}'
    print(status_msg.format(datetime.now().strftime('%b %d %Y @ %H:%M:%S'), 
                            url_info['url'], slice_method, loc_cmd))

    # Localize
    ecode = os.system(loc_cmd)
    if ecode != 0:
        err = 'LOCALIZATION ERROR: localization returned non-zero exit status for {}'
        exit(err.format(url_info['url']))

    # Index local slice (as sanity check)
    ecode = os.system(idx_cmd)
    if ecode != 0:
        err = 'INDEXING ERROR: index build failed for local slice {}'
        exit(err.format(url_info['basename_sliced']))


def _update_inputs_tsv(urls_dict, tsv_out):
    """
    Write output .tsv file with updated paths corresponding to local slices
    """

    if tsv_out in 'stdout -'.split():
        outfile = stdout
    else:
        outfile = open(tsv_out, 'w')

    for url_info in urls_dict.values():
        outfile.write('\t'.join([url_info['basename_sliced'], url_info['metadata']]) + '\n')

    outfile.close()


def slice_remote(urls_tsv, regions_bed, ref_fasta=None, local_suffix='local_slice', 
                 tsv_out='stdout'):
    """
    Localize slices of remote genomic data
    """

    # Parse input file
    urls_dict = _parse_remote_inputs_tsv(urls_tsv)

    # Make sure ref_fasta is provided if any CRAM URLs are supplied
    _check_cram_fasta_input(urls_dict, ref_fasta)

    # Authenticate GCP credentials (required for htslib streaming)
    _authenticate_remote_streaming()

    # Localize each slice
    for url_info in urls_dict.values():
        _localize_remote_slice(url_info, regions_bed, ref_fasta)

    # Write updated paths to tsv_out
    _update_inputs_tsv(urls_dict, tsv_out)

