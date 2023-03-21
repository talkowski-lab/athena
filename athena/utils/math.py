#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019- Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Mathematical utilities
"""


from scipy.stats import chi2, norm


def hwe_chisq(record):
    """
    Chi-squared test for deviation from Hardy-Weinberg equilibrium
    """

    # Get observed genotype counts
    aa_obs = record.info['N_HOMREF']
    Aa_obs = record.info['N_HET']
    AA_obs = record.info['N_HOMALT']
    total_N = aa_obs + Aa_obs + AA_obs

    # Get expected genotype counts
    a_freq = ((2 * aa_obs) + Aa_obs) / (2 * total_N)
    A_freq = ((2 * AA_obs) + Aa_obs) / (2 * total_N)
    aa_exp = (a_freq ** 2) * total_N
    Aa_exp = 2 * a_freq * A_freq * total_N
    AA_exp = (A_freq ** 2) * total_N

    # Calculate chi-square stats for each genotype
    aa_chisq = ((aa_obs - aa_exp) ** 2) / aa_exp
    Aa_chisq = ((Aa_obs - Aa_exp) ** 2) / Aa_exp
    AA_chisq = ((AA_obs - AA_exp) ** 2) / AA_exp
    chisq = aa_chisq + Aa_chisq + AA_chisq

    p = 1 - chi2.cdf(chisq, 1)

    return p


def ci_to_sd(ci_width, ci_pct):
        """
        Infer standard deviation from a confidence interval
        """

        # Note: ci_width measures the full width of the CI (lower to upper bound)
        # Thus, the total number of Z-scores covered in this interval is 2 * norm.ppf of the tail
        zscore_width = abs(2 * norm.ppf((1 - ci_pct) / 2))

        return ci_width / zscore_width

