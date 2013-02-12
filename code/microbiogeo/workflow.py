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
                               parse_anosim_permanova_results)
from microbiogeo.util import shuffle_dm, subsample_dm

def generate_distance_matrices(in_dir, out_dir, studies, metrics, num_shuffled,
                               num_subsets, tree_fp):
    # Process each depth in each study in parallel.
    c = Client()
    lview = c.load_balanced_view()
    lview.block = True

    create_dir(out_dir)

    per_study_depths = []
    for study in studies:
        for depth in studies[study]['depths']:
            per_study_depths.append((study, depth, metrics,
                                     studies[study]['categories'],
                                     studies[study]['group_sizes'],
                                     num_shuffled, num_subsets, in_dir,
                                     out_dir, tree_fp, shuffle_dm,
                                     subsample_dm))
    lview.map(generate_per_study_depth_dms, per_study_depths)

def main():
    in_dir = 'test_datasets'
    out_dir = 'test_grouping_analysis_output'
    tree_fp = join('test_datasets', 'overview', 'rep_set.tre')
    depth_descs = ['5_percent', '25_percent', '50_percent']
    studies = {
               'overview': {
                            'depths': [50, 100, 146],
                            'categories': ['Treatment'],
                            'group_sizes': [3, 4]
                           },
               'overview2': {
                             'depths': [50, 100, 146],
                             'categories': ['Treatment'],
                             'group_sizes': [3, 4]
                            }
              }
    metrics = ['euclidean', 'bray_curtis']
    grouping_methods = {
                        'adonis': parse_adonis_results,
                        'anosim': parse_anosim_permanova_results
                       }
    permutations = [99, 999]
    num_shuffled = 2
    num_subsets = 2

    generate_distance_matrices(in_dir, out_dir, studies, metrics, num_shuffled,
            num_subsets, tree_fp)


if __name__ == "__main__":
    main()
