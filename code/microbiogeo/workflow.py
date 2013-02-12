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

from glob import glob
from os.path import basename, exists, join, splitext

from IPython.parallel import Client

from qiime.util import create_dir

from microbiogeo.parallel import generate_per_study_depth_dms
from microbiogeo.parse import (parse_adonis_results,
                               parse_anosim_permanova_results,
                               parse_mantel_results,
                               parse_morans_i_results,
                               parse_partial_mantel_results)
from microbiogeo.util import run_command

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

def run_methods(in_dir, studies, grouping_methods, gradient_methods,
                permutations):
    """Runs each statistical method on each distance matrix."""
    # Process each compare_categories.py/compare_distance_matrices.py run in
    # parallel.
    c = Client()
    lview = c.load_balanced_view()
    lview.block = True

    jobs = []
    for study in studies:
        for depth in studies[study]['depths']:
            for method in grouping_methods:
                study_dir = join(in_dir, study)
                map_fp = join(study_dir, 'map.txt')
                depth_dir = join(study_dir, 'bdiv_even%d' % depth)
                dm_fps = glob(join(depth_dir, '*_dm*.txt'))

                for dm_fp in dm_fps:
                    for category in studies[study]['grouping_categories']:
                        for permutation in permutations:
                            dm_bn = basename(dm_fp)

                            if 'dm.txt' in dm_bn or 'dm_shuffled' in dm_bn or \
                               'dm_%s_gs' % category in dm_bn:
                                results_dir = join(depth_dir, '%s_%s_%s_%d' % (
                                    splitext(dm_bn)[0], method, category,
                                    permutation))

                                # Skip the job if the results dir exists. We'll
                                # assume it was created from a previous run.
                                if not exists(results_dir) or \
                                   len(listdir(results_dir)) == 0:
                                    cmd = ('compare_categories.py --method %s '
                                           '-n %d -i %s -m %s -c %s -o %s' % (
                                           method, permutation, dm_fp, map_fp,
                                           category, results_dir))
                                    jobs.append(cmd)

    lview.map(run_command, jobs)

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

    run_methods(out_dir, studies, grouping_methods, gradient_methods,
                permutations)


if __name__ == "__main__":
    main()
