#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Contains functions used in the ordination_correlation.py script."""

from qiime.format import format_p_value_for_num_iters
from qiime.parse import parse_coords, parse_mapping_file_to_dict
from qiime.pycogent_backports.test import correlation_test

CORRELATION_TYPES = ['pearson', 'spearman']

def compute_ordination_correlation(map_f, coord_f, category, axis=1,
                                   correlation_type='pearson',
                                   num_permutations=999):
    if correlation_type not in CORRELATION_TYPES:
        raise ValueError("Invalid correlation type '%s'. Must be one of %r." %
                         (correlation_type, CORRELATION_TYPES))
    if num_permutations < 0:
        raise ValueError("Invalid number of permutations: %d. Must be greater "
                         "than or equal to zero." % num_permutations)

    coords_samp_ids, coords, _, _ = parse_coords(coord_f)
    num_axes = len(coords[0])
    if axis < 1 or axis > num_axes:
        raise ValueError("Invalid axis number %d. Must be greater than zero "
                         "and less than or equal to the number of axes in the "
                         "input coordinates file (found %d axes)." %
                         (axis, num_axes))
    axis_data = coords[:, axis - 1]

    mdm, _ = parse_mapping_file_to_dict(map_f)
    gradient_data = []
    for samp_id in coords_samp_ids:
        if category not in mdm[samp_id]:
            raise ValueError("Category '%s' does not exist in the input "
                             "mapping file." % category)

        md_value = mdm[samp_id][category]
        try:
            md_value = float(md_value)
        except ValueError:
            raise ValueError("The category state '%s' could not be converted "
                             "to a number. All states in the '%s' category "
                             "must be numeric." % (md_value, category))
        gradient_data.append(md_value)

    corr_coeff, param_p_val, _, nonparam_p_val, _ = \
            correlation_test(axis_data, gradient_data, method=correlation_type,
                             permutations=num_permutations)

    if num_permutations > 0:
        nonparam_p_val = format_p_value_for_num_iters(nonparam_p_val,
                                                      num_permutations)
    else:
        nonparam_p_val = 'N/A'

    return corr_coeff, param_p_val, nonparam_p_val
