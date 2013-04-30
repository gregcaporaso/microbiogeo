#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Module for functionality meant to be run in parallel."""

from os.path import basename, join, splitext
from shutil import move

from qiime.filter import filter_samples_from_distance_matrix
from qiime.parse import parse_distmat
from qiime.util import create_dir

from microbiogeo.method import Mantel, MoransI
from microbiogeo.util import (choose_gradient_subsets, has_results,
                              run_command, shuffle_dm, subset_dm,
                              subset_groups)

# The parameters need to be wrapped in parens in order to work with map. We put
# these functions in their own module to avoid import issues when using
# IPython.parallel. See http://stackoverflow.com/a/12307741 for more details.
def generate_per_study_depth_dms((in_dir, out_dir, study, depth, metrics,
                                  grouping_categories, gradient_categories,
                                  group_sizes, subset_sizes, num_shuffled,
                                  num_subsets, tree_fp)):
    in_study_dir = join(in_dir, study)
    out_study_dir = join(out_dir, study)
    map_fp = join(out_study_dir, 'map.txt')

    full_otu_fp = join(in_study_dir, 'otu_table.biom')
    even_otu_fp = join(out_study_dir, 'otu_table_even%d.biom' % depth)
    bdiv_out_dir = join(out_study_dir, 'bdiv_even%d' % depth)

    metrics_param = ','.join(metrics)

    cmd = 'single_rarefaction.py -i %s -o %s -d %d' % (full_otu_fp,
                                                       even_otu_fp,
                                                       depth)
    run_command(cmd)

    cmd = 'beta_diversity.py -i %s -o %s -m %s -t %s' % (even_otu_fp,
                                                         bdiv_out_dir,
                                                         metrics_param,
                                                         tree_fp)
    run_command(cmd)

    # Rename each file to match QIIME's standard naming conventions of
    # distance matrices. Generate shuffled versions of each distance matrix.
    for metric in metrics:
        dm_fp = join(bdiv_out_dir, '%s_otu_table_even%d.txt' % (metric, depth))
        renamed_dm_fp = join(bdiv_out_dir, '%s_dm.txt' % metric)
        move(dm_fp, renamed_dm_fp)

        for i in range(1, num_shuffled + 1):
            renamed_dm_f = open(renamed_dm_fp, 'U')
            shuffled_dm_fp = join(bdiv_out_dir,
                                  '%s_dm_shuffled%d.txt' % (metric, i))
            shuffled_dm_f = open(shuffled_dm_fp, 'w')
            shuffled_dm_f.write(shuffle_dm(renamed_dm_f))
            shuffled_dm_f.close()
            renamed_dm_f.close()

        # Create distance matrices with subsets of groups.
        for category in grouping_categories:
            for group_size in group_sizes:
                for i in range(1, num_subsets + 1):
                    subset_dm_fp = join(bdiv_out_dir,
                            '%s_dm_%s_gs%d_%d.txt' % (metric, category,
                                                      group_size, i))
                    subset_dm_f = open(subset_dm_fp, 'w')
                    subset_dm_f.write(subset_groups(
                            open(renamed_dm_fp, 'U'), open(map_fp, 'U'),
                            category, group_size))
                    subset_dm_f.close()

        # Filter the keyboard study distance matrix to include only samples
        # taken from keys of subjects M2, M3, and M9 before creating subsets.
        # We only want to include these samples (not human subject fingertips)
        # because we want to see if keys that are closer to each other are
        # correlated with community similarity.
        if study == 'keyboard':
            filtered_dm_fp = join(bdiv_out_dir,
                    add_filename_suffix(renamed_dm_fp, '_keys_only'))
            cmd = ('filter_distance_matrix.py -i %s, -o %s -m %s -s '
                   '\'COMMON_NAME:keyboard;HOST_SUBJECT_ID:M2,M3,M9\'' % (
                       renamed_dm_fp, filtered_dm_fp, map_fp))
            run_command(cmd)
            renamed_dm_fp = filtered_dm_fp

        # Create subsets of each non-shuffled distance matrix.
        for subset_size in subset_sizes:
            for i in range(1, num_subsets + 1):
                subset_dm_fp = join(bdiv_out_dir,
                        '%s_dm_n%d_%d.txt' % (metric, subset_size, i))
                subset_dm_f = open(subset_dm_fp, 'w')
                subset_dm_f.write(
                        subset_dm(open(renamed_dm_fp, 'U'), subset_size))
                subset_dm_f.close()

# The following functions aren't actually run in parallel, but they are used to
# build commands that are run in parallel.
def build_grouping_method_commands(depth_dir, dm_fp, map_fp, method, category,
                                   permutations):
    cmds = []
    dm_bn = basename(dm_fp)

    for permutation in permutations:
        if 'dm.txt' in dm_bn or 'dm_shuffled' in dm_bn or \
           'dm_%s_gs' % category in dm_bn:
            results_dir = join(depth_dir, '%s_%s_%s_%d' % (splitext(dm_bn)[0],
                                                           method, category,
                                                           permutation))

            # Skip the job if the results dir exists. We'll assume it was
            # created from a previous run.
            if not has_results(results_dir):
                cmd = ('compare_categories.py --method %s -n %d -i %s -m %s '
                       '-c %s -o %s' % (method, permutation, dm_fp, map_fp,
                                        category, results_dir))
                cmds.append(cmd)
    return cmds

def build_gradient_method_commands(study_dir, depth_dir, dm_fp, map_fp, method,
                                   category, permutations):
    # Not all methods can be tested on the same variables/inputs. For example,
    # Moran's I cannot be tested on the keyboard study's key distances because
    # it is univariate, and partial Mantel can only be tested on the keyboard
    # study because that is the only study we have a control matrix for.
    cmds = []

    if method == 'mantel' or method == 'mantel_corr':
        for permutation in permutations:
            in_dm_fps = '%s,%s' % (dm_fp, join(study_dir,
                                               '%s_dm.txt' % category))
            results_dir = join(depth_dir,
                               '%s_%s_%s_%d' % (splitext(basename(dm_fp))[0],
                                                method, category, permutation))

            if not has_results(results_dir):
                cmd = ('compare_distance_matrices.py --method %s -n %d -i %s '
                       '-o %s' % (method, permutation, in_dm_fps, results_dir))
                cmds.append(cmd)
    elif method == 'morans_i':
        # Moran's I does not accept a number of permutations.
        results_dir = join(depth_dir,
                           '%s_%s_%s' % (splitext(basename(dm_fp))[0], method,
                                         category))

        if not has_results(results_dir):
            cmd = ('compare_categories.py --method %s -i %s -m %s -c %s '
                   '-o %s' % (method, dm_fp, map_fp, category, results_dir))
            cmds.append(cmd)

    return cmds

def build_gradient_method_keyboard_commands(study_dir, depth_dir, dm_fp,
                                            method, permutations):
    cmds = []

    if method in ('mantel', 'mantel_corr', 'partial_mantel'):
        for permutation in permutations:
            in_dm_fps = '%s,%s' % (dm_fp,
                                   join(study_dir,
                                        'euclidean_key_distances_dm.txt'))

            results_dir = join(depth_dir,
                               '%s_%s_%s_%d' % (splitext(basename(dm_fp))[0],
                                                method,
                                                'key_distance',
                                                permutation))

            if not has_results(results_dir):
                cmd = ('compare_distance_matrices.py --method %s -n %d -i %s '
                       '-o %s' % (method, permutation, in_dm_fps, results_dir))

                if method == 'partial_mantel':
                    control_dm_fp = join(study_dir,
                            'median_unifrac_individual_distances_dm.txt')
                    cmd += ' -c %s' % control_dm_fp

                cmds.append(cmd)

    return cmds

def build_best_method_commands(depth_dir, dm_fp, map_fp, env_vars):
    cmds = []

    results_dir = join(depth_dir,
                       '%s_%s' % (splitext(basename(dm_fp))[0], 'best'))

    if not has_results(results_dir):
        cmd = ('compare_categories.py --method best -i %s -m %s -c %s '
               '-o %s' % (dm_fp, map_fp, ','.join(env_vars), results_dir))
        cmds.append(cmd)

    return cmds

def build_grouping_method_sample_size_testing_commands(study_dir, dm_fp,
                                                       map_fp, metric,
                                                       category, methods,
                                                       subset_sizes,
                                                       num_subsets,
                                                       permutation):
    cmds = []

    for subset_size in subset_sizes:
        for subset_num in range(1, num_subsets + 1):
            map_f = open(map_fp, 'U')
            dm_f = open(dm_fp, 'U')
            subset_dm_fp = join(study_dir, '%s_dm_%s_gs%d_%d.txt' % (metric,
                    category, subset_size, subset_num))
            subset_dm_f = open(subset_dm_fp, 'w')
            subset_dm_f.write(
                    subset_groups(dm_f, map_f, category, subset_size))
            subset_dm_f.close()
            dm_f.close()
            map_f.close()

            for method in methods:
                results_dir = join(study_dir, '%s_dm_%s_gs%d_%d_%s' % (metric,
                        category, subset_size, subset_num, method.Name))

                if not has_results(results_dir):
                    cmd = ('compare_categories.py --method %s -n %d -i %s '
                           '-m %s -c %s -o %s' % (method.Name, permutation,
                                                  subset_dm_fp, map_fp,
                                                  category, results_dir))
                    cmds.append(cmd)

    return cmds

def build_gradient_method_sample_size_testing_commands(study_dir, dm_fp,
                                                       env_dm_fp, map_fp,
                                                       metric, category,
                                                       methods, subset_sizes,
                                                       num_subsets,
                                                       permutation):
    cmds = []

    dm_f = open(dm_fp, 'U')
    map_f = open(map_fp, 'U')
    sid_subsets = choose_gradient_subsets(dm_f, map_f, category, subset_sizes,
                                          num_subsets)
    map_f.close()

    dm_f.seek(0)
    dm = parse_distmat(dm_f)
    dm_f.close()

    results_idx = 0
    for subset_size in subset_sizes:
        for subset_num in range(1, num_subsets + 1):
            # Write filtered distance matrix subset.
            subset_dm_fp = join(study_dir, '%s_dm_%s_n%d_%d.txt' % (metric,
                    category, subset_size, subset_num))

            subset_dm_f = open(subset_dm_fp, 'w')
            subset_dm_f.write(filter_samples_from_distance_matrix(dm,
                              sid_subsets[results_idx], negate=True))
            subset_dm_f.close()

            for method in methods:
                results_dir = join(study_dir, '%s_dm_%s_n%d_%d_%s' % (metric,
                        category, subset_size, subset_num, method.Name))
                    
                if type(method) is Mantel:
                    in_dm_fps = '%s,%s' % (subset_dm_fp, env_dm_fp)

                    if not has_results(results_dir):
                        cmd = ('compare_distance_matrices.py --method %s '
                               '-n %d -i %s -o %s' % (method.Name, permutation,
                                                      in_dm_fps, results_dir))
                        cmds.append(cmd)
                elif type(method) is MoransI:
                    if not has_results(results_dir):
                        cmd = ('compare_categories.py --method %s -i %s -m %s '
                               '-c %s -o %s' % (method.Name, subset_dm_fp,
                                                map_fp, category, results_dir))
                        cmds.append(cmd)

            results_idx += 1

    return cmds
