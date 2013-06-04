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

from collections import defaultdict
from os import listdir
from os.path import basename, exists, join, splitext
from random import randint, sample

from biom.parse import parse_biom_table

from numpy import ceil, inf, mean, std

from qiime.filter import (filter_mapping_file_from_mapping_f,
                          filter_samples_from_otu_table)
from qiime.parse import (parse_mapping_file_to_dict, parse_mapping_file,
                         parse_coords, group_by_field)
from qiime.util import add_filename_suffix, create_dir, MetadataMap

from microbiogeo.method import (AbstractStatMethod, Adonis, Anosim, Best,
                                Dbrda, Mantel, MantelCorrelogram, MoransI,
                                Mrpp, PartialMantel, Permanova, Permdisp,
                                QiimeStatMethod, UnparsableFileError,
                                UnparsableLineError)
from microbiogeo.simulate import create_simulated_data_plots
from microbiogeo.util import (get_color_pool,
                              get_num_samples_in_distance_matrix,
                              get_num_samples_in_map, get_num_samples_in_table,
                              get_panel_label, get_simsam_rep_num, has_results,
                              run_command, run_parallel_jobs)

def generate_data(analysis_type, in_dir, out_dir, workflow, tree_fp):
    """Generates real and simulated data for each study.

    Distance matrices will be created at each even sampling depth and metric
    using the provided tree if necessary. Shuffled versions of each distance
    matrix will also be created, which can be used as negative controls.
    Additionally, simulated gradient or cluster data will be created at varying
    sample sizes and dissimilarity levels (using simsam.py).

    data_type should be either 'gradient' or 'cluster'.

    Will create the following (heavily nested) output directory structure:

    out_dir/
        study/
            depth/
                even depth otu table (.biom)
                real/
                    metric/
                        original/
                            map.txt
                            dm.txt
                            pc.txt
                            <category>_dm.txt (if gradient)
                        shuff_num
                            map.txt
                            dm.txt
                            pc.txt
                            <category>_dm.txt (if gradient)
                simulated/
                    category/
                        trial_num/
                            samp_size/
                                (optional) subset files dependent on samp_size
                                dissim/
                                    subset files/dirs dependent on samp_size
                                    metric/
                                        map.txt
                                        dm.txt
                                        pc.txt
                                        <category>_dm.txt (if gradient)
    """
    create_dir(out_dir)

    cmds = []
    for study in workflow:
        study_dir = join(out_dir, study)
        create_dir(study_dir)

        otu_table_fp = join(in_dir, study, 'otu_table.biom')
        map_fp = join(in_dir, study, 'map.txt')
        map_f = open(map_fp, 'U')

        for depth in workflow[study]['depths']:
            depth_dir = join(study_dir, '%d' % depth[0])
            create_dir(depth_dir)

            # Rarefy the table first since simsam.py's output tables will still
            # have even sampling depth and we don't want to lose simulated
            # samples after the fact.
            even_otu_table_fp = join(depth_dir, basename(otu_table_fp))

            if not exists(even_otu_table_fp):
                run_command('single_rarefaction.py -i %s -o %s -d %d;' % (
                        otu_table_fp, even_otu_table_fp, depth[0]))

            cmds.extend(_build_real_data_commands(analysis_type, depth_dir,
                    even_otu_table_fp, map_fp, tree_fp, workflow[study]))
            cmds.extend(_build_simulated_data_commands(analysis_type,
                    depth_dir, even_otu_table_fp, map_fp, tree_fp,
                    workflow[study]))

    run_parallel_jobs(cmds, run_command)

def _build_real_data_commands(analysis_type, out_dir, even_otu_table_fp,
                              map_fp, tree_fp, workflow):
    cmds = []

    data_type_dir = join(out_dir, 'real')
    create_dir(data_type_dir)

    for metric in workflow['metrics']:
        metric_dir = join(data_type_dir, metric[0])
        create_dir(metric_dir)

        orig_dir = join(metric_dir, 'original')
        create_dir(orig_dir)

        required_files = ['dm.txt', 'map.txt', 'pc.txt']
        if analysis_type == 'gradient':
            for category in workflow['categories']:
                required_files.append('%s_dm.txt' % category[0])

        has_orig_files = has_results(orig_dir, required_files=required_files)

        has_shuff_files = True
        for shuff_num in range(workflow['num_shuffled_trials']):
            shuff_num_dir = join(metric_dir, '%d' % shuff_num)
            has_shuff_files = has_results(shuff_num_dir, required_files)
            if not has_shuff_files:
                break

        if not (has_orig_files and has_shuff_files):
            cmds.append(_build_per_metric_real_data_commands(analysis_type,
                    metric_dir, even_otu_table_fp, map_fp, tree_fp, metric,
                    workflow['categories'], workflow['num_shuffled_trials']))
    return cmds

def _build_per_metric_real_data_commands(analysis_type, out_dir,
                                         even_otu_table_fp, map_fp, tree_fp,
                                         metric, categories,
                                         num_shuffled_trials):
    orig_dir = join(out_dir, 'original')

    cmd = ['beta_diversity.py -i %s -o %s -m %s -t %s' % (even_otu_table_fp, orig_dir, metric[0], tree_fp)]
    cmd.append('mv %s %s' % (join(orig_dir, '%s_%s.txt' % (metric[0], splitext(basename(even_otu_table_fp))[0])), join(orig_dir, 'dm.txt')))
    cmd.append('cp %s %s' % (map_fp, join(orig_dir, 'map.txt')))
    cmd.append('principal_coordinates.py -i %s -o %s' % (join(orig_dir, 'dm.txt'), join(orig_dir, 'pc.txt')))

    if analysis_type == 'gradient':
        for category in categories:
            cmd.append('distance_matrix_from_mapping.py -i %s -c %s -o %s' % (join(orig_dir, 'map.txt'), category[0], join(orig_dir, '%s_dm.txt' % category[0])))

    for shuff_num in range(num_shuffled_trials):
        shuff_num_dir = join(out_dir, '%d' % shuff_num)

        cmd.append('mkdir -p %s' % shuff_num_dir)
        cmd.append('shuffle_distance_matrix.py -i %s -o %s' % (join(orig_dir, 'dm.txt'), join(shuff_num_dir, 'dm.txt')))
        cmd.append('cp %s %s' % (join(orig_dir, 'map.txt'), join(shuff_num_dir, 'map.txt')))
        cmd.append('principal_coordinates.py -i %s -o %s' % (join(shuff_num_dir, 'dm.txt'), join(shuff_num_dir, 'pc.txt')))

        if analysis_type == 'gradient':
            for category in categories:
                cmd.append('distance_matrix_from_mapping.py -i %s -c %s -o %s' % (join(shuff_num_dir, 'map.txt'), category[0], join(shuff_num_dir, '%s_dm.txt' % category[0])))
    return ' && '.join(cmd)

def _build_simulated_data_commands(analysis_type, out_dir, even_otu_table_fp,
                                   map_fp, tree_fp, workflow):
    cmds = []

    data_type_dir = join(out_dir, 'simulated')
    create_dir(data_type_dir)

    num_samps = get_num_samples_in_table(even_otu_table_fp)

    for category in workflow['categories']:
        category_dir = join(data_type_dir, category[0])
        create_dir(category_dir)

        for trial_num in range(workflow['num_sim_data_trials']):
            trial_num_dir = join(category_dir, '%d' % trial_num)
            create_dir(trial_num_dir)

            for samp_size in workflow['sample_sizes']:
                samp_size_dir = join(trial_num_dir, '%d' % samp_size)
                create_dir(samp_size_dir)

                # Lots of duplicate code between these two blocks...
                # need to refactor and test.
                if samp_size <= num_samps:
                    simsam_rep_num = 1

                    subset_otu_table_fp = join(samp_size_dir, basename(even_otu_table_fp))
                    subset_map_fp = join(samp_size_dir, basename(map_fp))

                    if not has_results(samp_size_dir, required_files=[basename(subset_otu_table_fp), basename(subset_map_fp)]):
                        run_command('choose_data_subset.py -t %s -i %s -m %s -c %s -n %d -o %s' % (analysis_type, even_otu_table_fp, map_fp, category[0], samp_size, samp_size_dir))
                    assert get_num_samples_in_table(subset_otu_table_fp) == samp_size
                    assert get_num_samples_in_map(subset_map_fp) == samp_size

                    for d in workflow['dissim']:
                        dissim_dir = join(samp_size_dir, repr(d))
                        create_dir(dissim_dir)

                        simsam_map_fp = join(dissim_dir, add_filename_suffix(subset_map_fp, '_n%d_d%r' % (simsam_rep_num, d)))
                        simsam_otu_table_fp = join(dissim_dir, add_filename_suffix(subset_otu_table_fp, '_n%d_d%r' % (simsam_rep_num, d)))

                        # Check for simulated table/map and various
                        # distance matrices / coordinates files.
                        required_simsam_files = [basename(simsam_map_fp), basename(simsam_otu_table_fp)]
                        has_simsam_files = has_results(dissim_dir, required_files=required_simsam_files)

                        has_metric_files = True
                        for metric in workflow['metrics']:
                            required_metric_files = ['dm.txt', 'map.txt', 'pc.txt']
                            if analysis_type == 'gradient':
                                required_metric_files.append('%s_dm.txt' % category[0])

                            metric_dir = join(dissim_dir, metric[0])
                            has_metric_files = has_results(metric_dir, required_metric_files)
                            if not has_metric_files:
                                break

                        if not (has_simsam_files and has_metric_files):
                            cmd = ['simsam.py -i %s -t %s -o %s -d %r -n %d -m %s' % (subset_otu_table_fp, tree_fp, dissim_dir, d, simsam_rep_num, subset_map_fp)]

                            for metric in workflow['metrics']:
                                metric_dir = join(dissim_dir, metric[0])
                                create_dir(metric_dir)

                                if analysis_type == 'gradient':
                                    cmd.append('distance_matrix_from_mapping.py -i %s -c %s -o %s' % (simsam_map_fp, category[0], join(metric_dir, '%s_dm.txt' % category[0])))

                                cmd.append('beta_diversity.py -i %s -o %s -m %s -t %s' % (simsam_otu_table_fp, metric_dir, metric[0], tree_fp))
                                cmd.append('mv %s %s' % (join(metric_dir, '%s_%s.txt' % (metric[0], splitext(basename(simsam_otu_table_fp))[0])), join(metric_dir, 'dm.txt')))
                                cmd.append('cp %s %s' % (simsam_map_fp, join(metric_dir, 'map.txt')))
                                cmd.append('principal_coordinates.py -i %s -o %s' % (join(metric_dir, 'dm.txt'), join(metric_dir, 'pc.txt')))
                            cmds.append(' && '.join(cmd))
                else:
                    # We need to simulate more samples than we originally have.
                    simsam_rep_num = get_simsam_rep_num(samp_size, num_samps)

                    for d in workflow['dissim']:
                        dissim_dir = join(samp_size_dir, repr(d))
                        create_dir(dissim_dir)

                        simsam_map_fp = join(dissim_dir, add_filename_suffix(map_fp, '_n%d_d%r' % (simsam_rep_num, d)))
                        simsam_otu_table_fp = join(dissim_dir, add_filename_suffix(even_otu_table_fp, '_n%d_d%r' % (simsam_rep_num, d)))

                        required_simsam_files = [basename(simsam_map_fp), basename(simsam_otu_table_fp)]
                        has_simsam_files = has_results(dissim_dir, required_files=required_simsam_files)

                        required_subset_files = [basename(simsam_map_fp), basename(simsam_otu_table_fp)]
                        has_subset_files = has_results(join(dissim_dir, 'subset'), required_files=required_subset_files)

                        has_metric_files = True
                        for metric in workflow['metrics']:
                            required_metric_files = ['dm.txt', 'map.txt', 'pc.txt']
                            if analysis_type == 'gradient':
                                required_metric_files.append('%s_dm.txt' % category[0])

                            metric_dir = join(dissim_dir, metric[0])
                            has_metric_files = has_results(metric_dir, required_metric_files)
                            if not has_metric_files:
                                break

                        if not (has_simsam_files and has_subset_files and has_metric_files):
                            cmd = ['simsam.py -i %s -t %s -o %s -d %r -n %d -m %s' % (even_otu_table_fp, tree_fp, dissim_dir, d, simsam_rep_num, map_fp)]

                            subset_dir = join(dissim_dir, 'subset')
                            cmd.append('choose_data_subset.py -t %s -i %s -m %s -c %s -n %d -o %s' % (analysis_type, simsam_otu_table_fp, simsam_map_fp, category[0], samp_size, subset_dir))
                            subset_otu_table_fp = join(subset_dir, basename(simsam_otu_table_fp))
                            subset_map_fp = join(subset_dir, basename(simsam_map_fp))

                            for metric in workflow['metrics']:
                                metric_dir = join(dissim_dir, metric[0])
                                create_dir(metric_dir)

                                if analysis_type == 'gradient':
                                    cmd.append('distance_matrix_from_mapping.py -i %s -c %s -o %s' % (subset_map_fp, category[0], join(metric_dir, '%s_dm.txt' % category[0])))

                                cmd.append('beta_diversity.py -i %s -o %s -m %s -t %s' % (subset_otu_table_fp, metric_dir, metric[0], tree_fp))
                                cmd.append('mv %s %s' % (join(metric_dir, '%s_%s.txt' % (metric[0], splitext(basename(subset_otu_table_fp))[0])), join(metric_dir, 'dm.txt')))
                                cmd.append('cp %s %s' % (subset_map_fp, join(metric_dir, 'map.txt')))
                                cmd.append('principal_coordinates.py -i %s -o %s' % (join(metric_dir, 'dm.txt'), join(metric_dir, 'pc.txt')))
                            cmds.append(' && '.join(cmd))
    return cmds

def process_data(in_dir, workflow):
    """Run statistical methods over generated data.

    For real data, creates category and method dirs for original and shuffled
    data. Under each method dir, permutation dirs will also be
    created, e.g.:

    in_dir/
        ...
            category/
                method/
                    num_perms/
                        <method>_results.txt

    For simulated data, creates method dirs under metric dirs in in_dir, e.g.:

    in_dir/
        ...
            metric/
                method/
                    <method>_results.txt
    """
    # Process each compare_categories.py/compare_distance_matrices.py run in
    # parallel.
    cmds = []
    for study in workflow:
        study_dir = join(in_dir, study)

        for depth in workflow[study]['depths']:
            depth_dir = join(study_dir, '%d' % depth[0])

            cmds.extend(_build_real_data_methods_commands(depth_dir,
                    workflow[study]))
            cmds.extend(_build_simulated_data_methods_commands(depth_dir,
                    workflow[study]))

    run_parallel_jobs(cmds, run_command)

def _build_real_data_methods_commands(out_dir, workflow):
    cmds = []

    data_type_dir = join(out_dir, 'real')

    num_shuffled_trials = workflow['num_shuffled_trials']
    num_perms = workflow['num_real_data_perms']

    for metric in workflow['metrics']:
        metric_dir = join(data_type_dir, metric[0])

        dirs_to_process = ['original'] + map(str, range(num_shuffled_trials))
        for dir_to_process in dirs_to_process:
            dir_to_process = join(metric_dir, dir_to_process)

            dm_fp = join(dir_to_process, 'dm.txt')
            map_fp = join(dir_to_process, 'map.txt')

            for category in workflow['categories']:
                category_dir = join(dir_to_process, category[0])
                create_dir(category_dir)

                grad_dm_fp = join(dir_to_process, '%s_dm.txt' % category[0])

                for method in workflow['methods']:
                    if type(method) is Best or type(method) is PartialMantel:
                        continue

                    method_dir = join(category_dir, method.Name)
                    create_dir(method_dir)

                    if type(method) is MoransI:
                        if not has_results(method_dir):
                            cmds.append('compare_categories.py --method %s -i %s -m %s -c %s -o %s' % (method.Name, dm_fp, map_fp, category[0], method_dir))
                    else:
                        for perms in num_perms:
                            perms_dir = join(method_dir, '%d' % perms)
                            create_dir(perms_dir)

                            if not has_results(perms_dir):
                                if type(method) is Mantel or type(method) is MantelCorrelogram:
                                    in_dm_fps = ','.join((dm_fp, grad_dm_fp))
                                    cmds.append('compare_distance_matrices.py --method %s -n %d -i %s -o %s' % (method.Name, perms, in_dm_fps, perms_dir))
                                else:
                                    cmds.append('compare_categories.py --method %s -i %s -m %s -c %s -o %s -n %d' % (method.Name, dm_fp, map_fp, category[0], perms_dir, perms))

            if Best() in workflow['methods']:
                best_dir = join(dir_to_process, Best().Name)

                if not has_results(best_dir):
                    env_vars = ','.join(workflow['best_method_env_vars'])
                    cmds.append('compare_categories.py --method %s -i %s -m %s -c %s -o %s' % (Best().Name, dm_fp, map_fp, env_vars, best_dir))
    return cmds

def _build_simulated_data_methods_commands(out_dir, workflow):
    cmds = []

    data_type_dir = join(out_dir, 'simulated')

    num_sim_data_trials = workflow['num_sim_data_trials']
    num_sim_data_perms = workflow['num_sim_data_perms']

    for category in workflow['categories']:
        category_dir = join(data_type_dir, category[0])

        for trial_num in range(num_sim_data_trials):
            trial_num_dir = join(category_dir, '%d' % trial_num)

            for samp_size in workflow['sample_sizes']:
                samp_size_dir = join(trial_num_dir, '%d' % samp_size)

                for d in workflow['dissim']:
                    dissim_dir = join(samp_size_dir, repr(d))

                    for metric in workflow['metrics']:
                        metric_dir = join(dissim_dir, metric[0])

                        dm_fp = join(metric_dir, 'dm.txt')
                        map_fp = join(metric_dir, 'map.txt')
                        grad_dm_fp = join(metric_dir,
                                          '%s_dm.txt' % category[0])
                        assert get_num_samples_in_distance_matrix(dm_fp) == samp_size
                        assert get_num_samples_in_map(map_fp) == samp_size

                        for method in workflow['methods']:
                            if type(method) is Best or type(method) is PartialMantel:
                                continue
                            method_dir = join(metric_dir, method.Name)
                            create_dir(method_dir)

                            if not has_results(method_dir):
                                if type(method) is Mantel or type(method) is MantelCorrelogram:
                                    assert get_num_samples_in_distance_matrix(grad_dm_fp) == samp_size
                                    in_dm_fps = ','.join((dm_fp,
                                                          grad_dm_fp))
                                    cmds.append('compare_distance_matrices.py --method %s -n %d -i %s -o %s' % (method.Name, num_sim_data_perms, in_dm_fps, method_dir))
                                else:
                                    cmds.append('compare_categories.py --method %s -i %s -m %s -c %s -o %s -n %d' % (method.Name, dm_fp, map_fp, category[0], method_dir, num_sim_data_perms))
    return cmds

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

def main():
    test = True

    if test:
        in_dir = 'test_datasets'
        out_dir = 'test_output'
        out_gradient_dir = join(out_dir, 'gradient')
        out_cluster_dir = join(out_dir, 'cluster')
        tree_fp = join('test_datasets', 'overview', 'rep_set.tre')

        gradient_workflow = {
            'overview': {
                'categories': [('Gradient', 'Gradient Category')],
                'best_method_env_vars': ['Gradient'],
                'depths': [(146, '5_percent'), (148, '25_percent')],
                'metrics': [('unweighted_unifrac', 'Unweighted UniFrac'),
                            ('weighted_unifrac', 'Weighted UniFrac')],
                'num_real_data_perms': [99, 999],
                'num_sim_data_perms': 999,
                'dissim': [0.0, 0.001, 0.01, 0.1, 1.0, 10.0],
                'pcoa_dissim': [0.0, 0.001, 1.0, 10.0],
                'sample_sizes': [3, 5, 13],
                'pcoa_sample_size': 13,
                'num_sim_data_trials': 3,
                'num_shuffled_trials': 2,
                'methods': [Best(), Mantel(), MantelCorrelogram(), MoransI()]
            }
        }

        cluster_workflow = {
            'overview': {
                'categories': [('Treatment', 'Treatment Category',
                                {'Control': 'Control', 'Fast': 'Fast'})],
                'depths': [(146, '5_percent'), (148, '25_percent')],
                'metrics': [('unweighted_unifrac', 'Unweighted UniFrac'),
                            ('weighted_unifrac', 'Weighted UniFrac')],
                'num_real_data_perms': [99, 999],
                'num_sim_data_perms': 999,
                'dissim': [0.0, 0.001, 0.01, 0.1, 1.0, 10.0],
                'pcoa_dissim': [0.0, 0.001, 1.0, 10.0],
                'sample_sizes': [3, 5, 13],
                'pcoa_sample_size': 13,
                'num_sim_data_trials': 3,
                'num_shuffled_trials': 2,
                'methods': [Adonis(), Anosim()]
            }
        }

        heatmap_methods = {
            'grouping': [Adonis(), Anosim()],
            'gradient': [Mantel(), MoransI()]
        }
    else:
        in_dir = '../data'
        out_dir = 'microbiogeo_output'
        out_gradient_dir = join(out_dir, 'gradient')
        out_cluster_dir = join(out_dir, 'cluster')
        tree_fp = join('gg_otus_4feb2011', 'trees', 'gg_97_otus_4feb2011.tre')

        gradient_workflow = {
            '88_soils': {
                'categories': [('PH', 'pH'), ('LATITUDE', 'Latitude')],
                'best_method_env_vars': ['TOT_ORG_CARB', 'SILT_CLAY',
                    'ELEVATION', 'SOIL_MOISTURE_DEFICIT', 'CARB_NITRO_RATIO',
                    'ANNUAL_SEASON_TEMP', 'ANNUAL_SEASON_PRECPT', 'PH',
                    'CMIN_RATE', 'LONGITUDE', 'LATITUDE'
                ],
                'depths': [(400, '5_percent'), (580, '25_percent'),
                           (660, '50_percent')
                ],
                'metrics': [('unweighted_unifrac', 'Unweighted UniFrac'),
                            ('weighted_unifrac', 'Weighted UniFrac'),
                            ('bray_curtis', 'Bray-Curtis'),
                            ('euclidean', 'Euclidean')
                ],
                'num_real_data_perms': [99, 999],
                'num_sim_data_perms': 999,
                # dissim must all be floats!
                'dissim': [0.0, 0.001, 0.01, 0.1, 0.4, 0.7, 1.0, 10.0, 40.0,
                           70.0, 100.0],
                'pcoa_dissim': [0.0, 0.001, 1.0, 100.0],
                # sample_sizes must all be ints!
                'sample_sizes': [5, 10, 20, 40, 60, 80, 100, 150, 200, 300],
                'pcoa_sample_size': 150,
                'num_sim_data_trials': 10,
                'num_shuffled_trials': 5,
                'methods': [Best(), Mantel(), MantelCorrelogram(), MoransI()]
            },

            'gn': {
                'categories': [('LAYER', 'Layer')],
                'best_method_env_vars': ['LAYER', 'START_DEPTH', 'END_DEPTH',
                    'FERREDOXINSANDASSOCPROTEINS', 'FLAGELLA',
                    'PHOTOSYNTHESISRELATEDPROTEINS', 'CHEMOTAXIS',
                    'SUGARDEGRADATIONPATHWAYS', 'METHYLTRANSFERASE',
                    'ARYLSULFATASEAANDRELENZYMES', 'CHAPERONES',
                    'CYANOBACTERIALPROTEINDUF820'
                ],
                'depths': [(1276, '5_percent'), (1495, '25_percent'),
                           (1779, '50_percent')],
                'metrics': [('unweighted_unifrac', 'Unweighted UniFrac'),
                            ('weighted_unifrac', 'Weighted UniFrac'),
                            ('bray_curtis', 'Bray-Curtis'),
                            ('euclidean', 'Euclidean')
                ],
                'num_real_data_perms': [99, 999],
                'num_sim_data_perms': 999,
                'dissim': [0.0, 0.001, 0.01, 0.1, 0.4, 0.7, 1.0, 10.0, 40.0,
                           70.0, 100.0],
                'pcoa_dissim': [0.0, 0.001, 1.0, 100.0],
                'sample_sizes': [5, 10, 20, 40, 60, 80, 100, 150, 200, 300],
                'pcoa_sample_size': 150,
                'num_sim_data_trials': 10,
                'num_shuffled_trials': 5,
                'methods': [Best(), Mantel(), MantelCorrelogram(), MoransI()]
            }
        }

        cluster_workflow = {
            'keyboard': {
                'categories': [('HOST_SUBJECT_ID', 'Subject',
                                {'M2': 'Subject 1',
                                 'M3': 'Subject 2',
                                 'M9': 'Subject 3'})
                ],
                'depths': [(390, '5_percent'), (780, '25_percent'),
                           (1015, '50_percent')],
                'metrics': [('unweighted_unifrac', 'Unweighted UniFrac'),
                            ('weighted_unifrac', 'Weighted UniFrac'),
                            ('bray_curtis', 'Bray-Curtis'),
                            ('euclidean', 'Euclidean')
                ],
                'num_real_data_perms': [99, 999],
                'num_sim_data_perms': 999,
                'dissim': [0.0, 0.001, 0.01, 0.1, 0.4, 0.7, 1.0, 10.0, 40.0,
                           70.0, 100.0],
                # Good for plotting
                #'dissim': [0.0, 0.001, 0.01, 0.1, 1.0, 100.0],
                'pcoa_dissim': [0.0, 0.001, 1.0, 100.0],
                'sample_sizes': [5, 10, 20, 40, 60, 80, 100, 150, 200, 300],
                'pcoa_sample_size': 150,
                'num_sim_data_trials': 10,
                'num_shuffled_trials': 5,
                'methods': [Adonis(), Anosim(), Mrpp(), Permanova(), Dbrda(),
                            Permdisp()]
            },

            'whole_body': {
                'categories': [('BODY_SITE_COARSE', 'Body Site',
                                {'gut': 'Gut',
                                 'oral': 'Oral cavity',
                                 'skin': 'Skin'}),
                               ('SEX', 'Sex',
                                {'female': 'Female', 'male': 'Male'})
                ],
                'depths': [(575, '5_percent'), (877, '25_percent'),
                           (1110, '50_percent')],
                'metrics': [('unweighted_unifrac', 'Unweighted UniFrac'),
                            ('weighted_unifrac', 'Weighted UniFrac'),
                            ('bray_curtis', 'Bray-Curtis'),
                            ('euclidean', 'Euclidean')
                ],
                'num_real_data_perms': [99, 999],
                'num_sim_data_perms': 999,
                'dissim': [0.0, 0.001, 0.01, 0.1, 0.4, 0.7, 1.0, 10.0, 40.0,
                           70.0, 100.0],
                'pcoa_dissim': [0.0, 0.001, 1.0, 100.0],
                'sample_sizes': [5, 20, 40, 80, 140, 220, 320, 420, 520, 600],
                'pcoa_sample_size': 150,
                'num_sim_data_trials': 10,
                'num_shuffled_trials': 5,
                'methods': [Adonis(), Anosim(), Mrpp(), Permanova(), Dbrda(),
                            Permdisp()]
            }
        }

        heatmap_methods = {
            'grouping': [Adonis(), Dbrda(), Mrpp(), Permanova(), Anosim()],
            'gradient': [Mantel(), MoransI()]
        }

    generate_data('gradient', in_dir, out_gradient_dir, gradient_workflow,
                  tree_fp)
    generate_data('cluster', in_dir, out_cluster_dir, cluster_workflow, tree_fp)

    process_data(out_gradient_dir, gradient_workflow)
    process_data(out_cluster_dir, cluster_workflow)

    #create_real_data_summary_tables('gradient', out_gradient_dir,
    #                                gradient_workflow)
    #create_real_data_summary_tables('cluster', out_cluster_dir,
    #                                cluster_workflow)

    #create_simulated_data_plots('gradient', out_gradient_dir,
    #                            gradient_workflow)
    #create_simulated_data_plots('cluster', out_cluster_dir, cluster_workflow)

    #create_method_heatmaps()

if __name__ == "__main__":
    main()
