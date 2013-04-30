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
from os.path import exists, join
from shutil import copy

from numpy import median

from qiime.util import create_dir

from microbiogeo.format import (create_method_comparison_heatmaps,
                                create_results_summary_tables,
                                create_sample_size_plots)
from microbiogeo.parallel import (build_best_method_commands,
        build_gradient_method_commands,
        build_gradient_method_keyboard_commands,
        build_gradient_method_sample_size_testing_commands,
        build_grouping_method_commands,
        build_grouping_method_sample_size_testing_commands,
        generate_per_study_depth_dms)
from microbiogeo.method import (AbstractStatMethod, Adonis, Anosim, Best,
                                Dbrda, Mantel, MantelCorrelogram, MoransI,
                                Mrpp, PartialMantel, Permanova, Permdisp,
                                QiimeStatMethod, UnparsableFileError,
                                UnparsableLineError)
from microbiogeo.util import run_command, run_parallel_jobs, StatsResults

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
    create_dir(out_dir)

    # Process each depth in each study in parallel.
    per_study_depths = []
    for study in studies:
        # Prep per-study output directories before running each depth in
        # parallel.
        in_study_dir = join(in_dir, study)
        out_study_dir = join(out_dir, study)
        create_dir(out_study_dir)
        copy(join(in_study_dir, 'map.txt'), out_study_dir)
        map_fp = join(out_study_dir, 'map.txt')

        # Create distance matrices from environmental variables in mapping
        # file. These are independent of sampling depth and metric, so we only
        # need to create them once for each study. Again, keyboard is unique in
        # that we cannot easily create a key distance matrix from the mapping
        # file. We'll use one that has been precalculated.
        for category in studies[study]['gradient_categories']:
            env_dm_fp = join(out_study_dir, '%s_dm.txt' % category)

            cmd = ('distance_matrix_from_mapping.py -i %s -c %s -o %s' % (
                   map_fp, category, env_dm_fp))
            run_command(cmd)

        if study == 'keyboard':
            key_dm_fp = join(in_study_dir, 'euclidean_key_distances_dm.txt')
            copy(key_dm_fp, out_study_dir)

            indiv_dm_fp = join(in_study_dir,
                               'median_unifrac_individual_distances_dm.txt')
            copy(indiv_dm_fp, out_study_dir)

        for depth in studies[study]['depths']:
            per_study_depths.append((in_dir, out_dir, study, depth, metrics,
                                     studies[study]['grouping_categories'],
                                     studies[study]['gradient_categories'],
                                     studies[study]['group_sizes'],
                                     studies[study]['subset_sizes'],
                                     num_shuffled, num_subsets, tree_fp))

    run_parallel_jobs(per_study_depths, generate_per_study_depth_dms)

def run_methods(in_dir, studies, methods, permutations):
    """Runs each statistical method on each distance matrix."""
    # Process each compare_categories.py/compare_distance_matrices.py run in
    # parallel.
    jobs = []

    for study in studies:
        best_method_env_vars = studies[study]['best_method_env_vars']

        for depth in studies[study]['depths']:
            for method_type in methods:
                for method in methods[method_type]:
                    study_dir = join(in_dir, study)
                    map_fp = join(study_dir, 'map.txt')
                    depth_dir = join(study_dir, 'bdiv_even%d' % depth)
                    dm_fps = glob(join(depth_dir, '*_dm*.txt'))

                    for dm_fp in dm_fps:
                        if method_type == 'grouping':
                            for category in \
                                    studies[study]['grouping_categories']:
                                jobs.extend(build_grouping_method_commands(
                                        depth_dir, dm_fp, map_fp, method.Name,
                                        category, permutations))
                        elif method_type == 'gradient':
                            for category in \
                                    studies[study]['gradient_categories']:
                                jobs.extend(build_gradient_method_commands(
                                        study_dir, depth_dir, dm_fp, map_fp,
                                        method.Name, category, permutations))

                            # Handle special cases here.
                            if study == 'keyboard':
                                jobs.extend(
                                    build_gradient_method_keyboard_commands(
                                            study_dir, depth_dir, dm_fp,
                                            method.Name, permutations))

                            if type(method) is Best and best_method_env_vars:
                                jobs.extend(
                                    build_best_method_commands(depth_dir,
                                        dm_fp, map_fp, best_method_env_vars))
                        else:
                            raise ValueError("Unknown method type '%s'." %
                                             method_type)

    run_parallel_jobs(jobs, run_command)

def summarize_results(in_dir, out_dir, studies, methods, heatmap_methods,
                      depth_descs, metrics, permutations, num_shuffled,
                      num_subsets):
    """Summarizes the results of the various method runs.

    Effect size statistics and p-values are collected for each of the tests
    that were run and summary tables are created, one for each sampling depth /
    metric combination (separate tables for grouping and gradient analysis
    methods). These tables are written to out_dir.
    """
    results = _collate_results(in_dir, studies, methods, depth_descs, metrics,
                               permutations, num_shuffled, num_subsets)

    create_results_summary_tables(results, out_dir)

    create_method_comparison_heatmaps(results, heatmap_methods, out_dir)

def _collate_results(in_dir, studies, methods, depth_descs, metrics,
                     permutations, num_shuffled, num_subsets):
    results = {}

    for depth_idx, depth_desc in enumerate(depth_descs):
        depth_res = {}

        for metric in metrics:
            metric_res = {}

            for method_type in methods:
                method_type_res = {}

                for method in methods[method_type]:
                    if type(method) in (MantelCorrelogram, Best):
                        # Completely ignore Mantel correlogram and BEST (for
                        # now at least). Mantel correlogram is hard to
                        # summarize because it produces a correlogram and many
                        # Mantel statistics for each distance class. We'll need
                        # to look at those results by hand and summarize them
                        # in the paper. The same holds for BEST: though it does
                        # not create a visual plot, it does not provide
                        # p-values. It mainly tells you which environmental
                        # variables best correlate with the community data.
                        continue

                    method_res = {}

                    for study in studies:
                        study_res = {}

                        # Figure out what our actual depth is for the study,
                        # and what subset sizes we used.
                        depth = studies[study]['depths'][depth_idx]

                        if method_type == 'grouping':
                            subset_sizes = studies[study]['group_sizes']
                            categories = studies[study]['grouping_categories']
                        elif method_type == 'gradient':
                            subset_sizes = studies[study]['subset_sizes']
                            categories = studies[study]['gradient_categories']

                            # Add our fictional 'key_distance' category, which
                            # isn't actually a category (i.e. not in a mapping
                            # file), but can be treated the same way as the
                            # others in this case.
                            if study == 'keyboard':
                                categories = categories + ['key_distance']
                        else:
                            raise ValueError("Unknown method type '%s'." %
                                             method_type)

                        for category in categories:
                            category_res = {}
                            full_results = StatsResults()
                            shuffled_results = StatsResults()
                            ss_results = [StatsResults()
                                          for i in range(len(subset_sizes))]

                            # Moran's I does not use permutations.
                            if type(method) is MoransI:
                                _collate_category_results(full_results,
                                        shuffled_results, ss_results, in_dir,
                                        study, depth, metric, method_type,
                                        method, category, subset_sizes,
                                        num_shuffled, num_subsets,
                                        permutation=None)
                            else:
                                for permutation in permutations:
                                    _collate_category_results(full_results,
                                            shuffled_results, ss_results,
                                            in_dir, study, depth, metric,
                                            method_type, method, category,
                                            subset_sizes, num_shuffled,
                                            num_subsets,
                                            permutation=permutation)

                            category_res['full'] = full_results
                            category_res['shuffled'] = shuffled_results
                            category_res['subsampled'] = ss_results

                            study_res[category] = category_res
                        method_res[study] = study_res
                    method_type_res[method.Name] = method_res
                metric_res[method_type] = method_type_res
            depth_res[metric] = metric_res
        results[depth_desc] = depth_res

    return results

def _collate_category_results(full_results, shuffled_results, ss_results,
                              in_dir, study, depth, metric, method_type,
                              method, category, subset_sizes, num_shuffled,
                              num_subsets, permutation=None):
    depth_dir = join(in_dir, study, 'bdiv_even%d' % depth)

    # Collect results for full distance matrices.
    results_dir = join(depth_dir, '%s_dm_%s_%s' % (metric, method.Name,
                                                   category))
    if permutation is not None:
        results_dir += '_%d' % permutation

    results_fp = join(results_dir, '%s_results.txt' % method.Name)

    # We will not always have results for every combination of parameters (e.g.
    # partial Mantel).
    if exists(results_fp):
        full_res_f = open(results_fp, 'U')
        full_es, full_p_val = method.parse(full_res_f)
        full_res_f.close()
        full_results.addResult(full_es, full_p_val)

    # Collect results for shuffled distance matrices.
    shuff_ess = []
    shuff_p_vals = []
    for shuff_num in range(1, num_shuffled + 1):
        results_dir = join(depth_dir, '%s_dm_shuffled%d_%s_%s' % (metric,
                                                                  shuff_num,
                                                                  method.Name,
                                                                  category))
        if permutation is not None:
            results_dir += '_%d' % permutation

        results_fp = join(results_dir, '%s_results.txt' % method.Name)

        if exists(results_fp):
            shuff_res_f = open(results_fp, 'U')
            shuff_es, shuff_p_val = method.parse(shuff_res_f)
            shuff_res_f.close()
            shuff_ess.append(shuff_es)
            shuff_p_vals.append(shuff_p_val)

    if shuff_ess and shuff_p_vals:
        shuffled_results.addResult(median(shuff_ess), median(shuff_p_vals))

    # Collect results for subset distance matrices.
    for subset_size_idx, subset_size in enumerate(subset_sizes):
        ss_ess = []
        ss_p_vals = []

        for ss_num in range(1, num_subsets + 1):
            if method_type == 'grouping':
                subset_path = '%s_dm_%s_gs%d_%d_%s_%s' % (metric, category,
                                                          subset_size, ss_num,
                                                          method.Name, category)
            elif method_type == 'gradient':
                subset_path = '%s_dm_n%d_%d_%s_%s' % (metric, subset_size,
                                                      ss_num, method.Name,
                                                      category)
            else:
                raise ValueError("Unknown method type '%s'." % method_type)

            results_dir = join(depth_dir, subset_path)
            if permutation is not None:
                results_dir += '_%d' % permutation

            results_fp = join(results_dir, '%s_results.txt' % method.Name)

            if exists(results_fp):
                ss_res_f = open(results_fp, 'U')
                ss_es, ss_p_val = method.parse(ss_res_f)
                ss_res_f.close()
                ss_ess.append(ss_es)
                ss_p_vals.append(ss_p_val)

        if ss_ess and ss_p_vals:
            ss_results[subset_size_idx].addResult(
                    median(ss_ess), median(ss_p_vals))

def run_sample_size_tests(in_dir, out_dir, sample_size_tests):
    """Tests the methods on subsets of the original sample groups.

    For grouping analysis methods:

    Samples will be chosen randomly without replacement to generate groups of
    samples at the specified subset size.

    For gradient analysis methods:

    Samples will be chosen along the gradient for each subset (which is what a
    researcher might do instead of randomly picking samples in order to test a
    gradient).

    Plots will be generated with subset size on the x-axis and test statistic
    on the y-axis. This will allow us to see if there is a cutoff/threshold
    based on the number of samples (i.e. where things start to stabilize in
    terms of the number of samples).
    """
    create_dir(out_dir)

    jobs = []
    for method_type in sample_size_tests:
        study = sample_size_tests[method_type]['study']
        depth = sample_size_tests[method_type]['depth']
        metric = sample_size_tests[method_type]['metric']
        categories = sample_size_tests[method_type]['categories']
        subset_sizes = sample_size_tests[method_type]['subset_sizes']
        num_subsets = sample_size_tests[method_type]['num_subsets']
        permutation = sample_size_tests[method_type]['permutation']
        methods = sample_size_tests[method_type]['methods']

        out_method_type_dir = join(out_dir, method_type)
        out_study_dir = join(out_method_type_dir, study)
        create_dir(out_method_type_dir)
        create_dir(out_study_dir)

        map_fp = join(in_dir, study, 'map.txt')
        dm_fp = join(in_dir, study, 'bdiv_even%d' % depth,
                     '%s_dm.txt' % metric)

        for category in categories:
            if method_type == 'grouping':
                jobs.extend(build_grouping_method_sample_size_testing_commands(
                        out_study_dir, dm_fp, map_fp, metric, category,
                        methods, subset_sizes, num_subsets, permutation))
            elif method_type == 'gradient':
                env_dm_fp = join(in_dir, study, '%s_dm.txt' % category)

                jobs.extend(build_gradient_method_sample_size_testing_commands(
                        out_study_dir, dm_fp, env_dm_fp, map_fp, metric,
                        category, methods, subset_sizes, num_subsets,
                        permutation))
            else:
                raise ValueError("Unknown method type '%s'." % method_type)

    run_parallel_jobs(jobs, run_command)

    create_sample_size_plots(out_dir, out_dir, sample_size_tests)

def main():
    test = True

    if test:
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
        methods = {
            'grouping': [Adonis(), Anosim()],
            'gradient': [Best(), Mantel(), MantelCorrelogram(), MoransI(),
                         PartialMantel()]
        }

        heatmap_methods = {
            'grouping': [Adonis(), Anosim()],
            'gradient': [Mantel(), MoransI()]
        }

        permutations = [99, 999]
        num_shuffled = 2
        num_subsets = 2

        # For sample size testing.
        sample_size_tests = {
            'grouping': {
                'study': 'overview',
                'depth': 100,
                'metric': 'bray_curtis',
                'subset_sizes': [3, 4],
                'num_subsets': 10,
                'permutation': 999,
                'categories': {
                    'Treatment': ['b', 'Treatment'],
                    'DOB': ['r', 'Date of birth']
                },
                'methods': [Adonis(), Anosim()]
            },

            'gradient': {
                'study': 'overview',
                'depth': 146,
                'metric': 'euclidean',
                'subset_sizes': [3, 4],
                'num_subsets': 20,
                'permutation': 999,
                'categories': {
                    'DOB': ['b', 'Date of birth']
                },
                'methods': [Mantel(), MoransI()]
            }
        }
    else:
        in_dir = 'datasets'
        out_dir = 'microbiogeo_output'
        tree_fp = join('gg_otus_4feb2011', 'trees', 'gg_97_otus_4feb2011.tre')
        depth_descs = ['5_percent', '25_percent', '50_percent']
        studies = {
                   '88_soils': {
                                'depths': [400, 580, 660],
                                'grouping_categories': ['ENV_BIOME'],
                                'gradient_categories': ['LATITUDE', 'PH'],
                                'group_sizes': [5, 10, 20, 40],
                                'subset_sizes': [10, 20, 30],
                                'best_method_env_vars': ['TOT_ORG_CARB',
                                    'SILT_CLAY', 'ELEVATION',
                                    'SOIL_MOISTURE_DEFICIT',
                                    'CARB_NITRO_RATIO', 'ANNUAL_SEASON_TEMP',
                                    'ANNUAL_SEASON_PRECPT', 'PH', 'CMIN_RATE',
                                    'LONGITUDE', 'LATITUDE']
                               }, 
                   'glen_canyon': {
                                   'depths': [15000, 29000, 53000],
                                   'grouping_categories': ['CurrentlyWet'],
                                   'gradient_categories': ['estimated_years_since_submerged_for_plotting'],
                                   'group_sizes': [5, 10, 20, 40],
                                   'subset_sizes': [10, 20, 30],
                                   'best_method_env_vars': ['sample_pH',
                                       'estimated_years_since_submerged_for_plotting',
                                       'Month', 'Day', 'Year',
                                       'days_since_epoch', 'Hour', 'Replicate',
                                       'DNA.I.D.No.']
                                  },
                   'keyboard': {
                                'depths': [390, 780, 1015],
                                'grouping_categories': ['HOST_SUBJECT_ID'],
                                'gradient_categories': [],
                                'group_sizes': [5, 10, 20, 40],
                                'subset_sizes': [10, 20, 30],
                                'best_method_env_vars': []
                               },
                   'whole_body': {
                                  'depths': [575, 877, 1110],
                                  'grouping_categories': ['BODY_SITE', 'SEX'],
                                  'gradient_categories': [],
                                  'group_sizes': [5, 10, 20, 40],
                                  'subset_sizes': [10, 20, 30],
                                  'best_method_env_vars': []
                                 }
                  }

        metrics = ['euclidean', 'bray_curtis', 'weighted_unifrac',
                   'unweighted_unifrac']

        methods = {
            'grouping': [Adonis(), Anosim(), Mrpp(), Permanova(), Dbrda(),
                         Permdisp()],
            'gradient': [Best(), Mantel(), MantelCorrelogram(), MoransI(),
                         PartialMantel()]
        }

        heatmap_methods = {
            'grouping': [Adonis(), Anosim(), Mrpp(), Permanova(), Dbrda()],
            'gradient': [Mantel(), MoransI()]
        }

        permutations = [99, 999]
        num_shuffled = 5
        num_subsets = 5

        # For sample size testing.
        sample_size_tests = {
            'grouping': {
                'study': 'whole_body',
                'depth': 575,
                'metric': 'unweighted_unifrac',
                'subset_sizes': [5, 10, 20, 40, 60, 80],
                'num_subsets': 10,
                'permutation': 999,
                'categories': {
                    'BODY_SITE': ['b', 'Body site'],
                    'SEX': ['r', 'Sex']
                },
                'methods': [Adonis(), Anosim(), Mrpp(), Permanova(), Dbrda()]
            },

            'gradient': {
                'study': '88_soils',
                'depth': 400,
                'metric': 'unweighted_unifrac',
                'subset_sizes': [5, 10, 20, 40, 60, 80],
                'num_subsets': 10,
                'permutation': 999,
                'categories': {
                    'PH': ['b', 'pH'],
                    'LATITUDE': ['r', 'Latitude']
                },
                'methods': [Mantel(), MoransI()]
            }
        }

    generate_distance_matrices(in_dir, out_dir, studies, metrics, num_shuffled,
            num_subsets, tree_fp)

    run_methods(out_dir, studies, methods, permutations)

    summarize_results(out_dir, out_dir, studies, methods, heatmap_methods,
                      depth_descs, metrics, permutations, num_shuffled,
                      num_subsets)

    run_sample_size_tests(out_dir, join(out_dir, 'sample_size_testing_output'),
                          sample_size_tests)


if __name__ == "__main__":
    main()
