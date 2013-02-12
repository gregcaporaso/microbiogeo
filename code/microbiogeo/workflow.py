#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Module for executing various workflows/analyses."""

from os.path import join

from IPython.parallel import Client

from qiime.util import create_dir

from microbiogeo.parallel import generate_per_study_depth_dms
from microbiogeo.parse import (parse_adonis_results,
                               parse_anosim_permanova_results,
                               parse_mantel_results,
                               parse_morans_i_results,
                               parse_partial_mantel_results)

def generate_distance_matrices(in_dir, out_dir, studies, metrics, num_shuffled,
                               num_subsets, tree_fp):
    """Generates distance matrices for each study.

    Distance matrices will be created at each even sampling depth and metric
    using the provided tree if necessary. Shuffled versions of each distance
    matrix will also be created, which can be used as negative controls.

    In addition, subsets of each distance matrix will be generated (based on a
    category's sample groupings or the entire distance matrix). These can be
    used to test how the methods perform on different study sizes.
    """
    # Process each depth in each study in parallel.
    c = Client()
    lview = c.load_balanced_view()
    lview.block = True

    create_dir(out_dir)

    per_study_depths = []
    for study in studies:
        for depth in studies[study]['depths']:
            per_study_depths.append((in_dir, out_dir, study, depth, metrics,
                                     studies[study]['grouping_categories'],
                                     studies[study]['gradient_categories'],
                                     studies[study]['group_sizes'],
                                     studies[study]['subset_sizes'],
                                     num_shuffled, num_subsets, tree_fp))
    lview.map(generate_per_study_depth_dms, per_study_depths)

def main():
    in_dir = 'test_datasets'
    out_dir = 'test_output'
    tree_fp = join('test_datasets', 'overview', 'rep_set.tre')
    depth_descs = ['5_percent', '25_percent', '50_percent']
    studies = {
               'overview': {
                            'depths': [50, 100, 146],
                            'grouping_categories': ['Treatment'],
                            'gradient_categories': ['DOB'],
                            'group_sizes': [3, 4],
                            'subset_sizes': [3, 4],
                            'best_method_env_vars': ['DOB']
                           },
               'overview2': {
                             'depths': [50, 100, 146],
                             'grouping_categories': ['Treatment'],
                             'gradient_categories': [],
                             'group_sizes': [3, 4],
                             'subset_sizes': [],
                             'best_method_env_vars': []
                            }
              }
    metrics = ['euclidean', 'bray_curtis']
    grouping_methods = {
            'adonis': parse_adonis_results,
            'anosim': parse_anosim_permanova_results
    }
    gradient_methods = {
            'mantel': parse_mantel_results,
            'mantel_corr': None,
            'morans_i': parse_morans_i_results,
            'partial_mantel': parse_partial_mantel_results
    }

    permutations = [99, 999]
    num_shuffled = 2
    num_subsets = 2

    generate_distance_matrices(in_dir, out_dir, studies, metrics, num_shuffled,
            num_subsets, tree_fp)


if __name__ == "__main__":
    main()
