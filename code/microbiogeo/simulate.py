#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Module with functionality for simulating gradient and cluster data."""

from biom.parse import parse_biom_table
from collections import defaultdict
from numpy import ceil, inf, mean, std
from os import listdir
from os.path import basename, exists, join, splitext
from matplotlib.lines import Line2D
from matplotlib.pyplot import (cm, figure, legend, figlegend, scatter, subplot,
                               title, xlim, xlabel, ylabel, xticks, yticks)
from qiime.filter import (filter_mapping_file_from_mapping_f,
                          filter_samples_from_otu_table)
from qiime.parse import (parse_mapping_file_to_dict, parse_mapping_file,
                         parse_coords, group_by_field)
from qiime.util import add_filename_suffix, create_dir, MetadataMap
from random import randint, sample

from microbiogeo.parse import (parse_adonis_results,
                               parse_anosim_permanova_results,
                               parse_dbrda_results,
                               parse_mantel_results,
                               parse_morans_i_results,
                               parse_mrpp_results)
from microbiogeo.util import (get_color_pool, get_num_samples, has_results,
                              run_command, run_parallel_jobs)

class InvalidSubsetSize(Exception):
    pass

def choose_cluster_subsets(otu_table_f, map_f, category, num_total_samples):
    otu_table = parse_biom_table(otu_table_f)
    metadata_map = MetadataMap.parseMetadataMap(map_f)

    # Dirty... :(
    try:
        map_f.seek(0)
    except AttributeError:
        pass

    if num_total_samples > len(otu_table.SampleIds):
        raise InvalidSubsetSize("Too many total samples (%d) were specified "
                                "as a subset size. There are only %d total "
                                "samples to choose a subset from." %
                                (num_total_samples, len(otu_table.SampleIds)))

    category_map = defaultdict(list)
    for samp_id in metadata_map.SampleIds:
        # Mapping files can have more samples than OTU tables.
        if samp_id in otu_table.SampleIds:
            category_val = metadata_map.getCategoryValue(samp_id, category)
            category_map[category_val].append(samp_id)

    samp_ids_to_keep, extra_samps = _choose_items_from_clusters(
            category_map, otu_table.SampleIds, num_total_samples)
    samp_ids_to_keep.extend(extra_samps)

    assert len(samp_ids_to_keep) == num_total_samples, \
           "%d != %d" % (len(samp_ids_to_keep), num_total_samples)
    assert len(samp_ids_to_keep) == len(set(samp_ids_to_keep)), \
           "Duplicate sample IDs in subset"

    return (filter_samples_from_otu_table(otu_table, samp_ids_to_keep, 0, inf),
            filter_mapping_file_from_mapping_f(map_f, samp_ids_to_keep))

def _choose_items_from_clusters(category_map, all_samp_ids, num_total_samples):
    # How many clusters we have.
    num_cat_states = len(category_map)

    # The number of samples to choose from each cluster.
    cluster_subset_size = num_total_samples // num_cat_states

    # Sort category states to facilitate unit testing.
    samp_ids_to_keep = []
    for category_val, samp_ids in sorted(category_map.items()):
        samp_ids_to_keep.extend(
                sample(samp_ids, min(cluster_subset_size, len(samp_ids))))

    remaining_samp_ids = set(all_samp_ids) - set(samp_ids_to_keep)

    # The number of remaining samples that we need to randomly choose from
    # (regardless of what cluster they are in) in order to meet our
    # num_total_samples quota.
    num_remaining_samps = num_total_samples - len(samp_ids_to_keep)

    return samp_ids_to_keep, sample(remaining_samp_ids, num_remaining_samps)

def choose_gradient_subset(otu_table_f, map_f, category, num_total_samples):
    otu_table = parse_biom_table(otu_table_f)
    mdm, _ = parse_mapping_file_to_dict(map_f)

    try:
        map_f.seek(0)
    except AttributeError:
        pass

    if num_total_samples > len(otu_table.SampleIds):
        raise InvalidSubsetSize("Too many total samples (%d) were specified "
                                "as a gradient subset size. There are only %d "
                                "total samples to choose a subset from." %
                                (num_total_samples, len(otu_table.SampleIds)))

    # Only keep the sample IDs that are in both the mapping file and OTU table.
    # Sort the samples according to the gradient category.
    samp_ids = [(samp_id, float(metadata[category]))
                for samp_id, metadata in mdm.items()
                if samp_id in otu_table.SampleIds]
    samp_ids.sort(key=lambda samp_id: samp_id[1])

    samp_ids_to_keep = [samp_id[0] for samp_id in
                        _choose_items_from_bins(samp_ids, num_total_samples)]

    assert len(samp_ids_to_keep) == num_total_samples, \
           "%d != %d" % (len(samp_ids_to_keep), num_total_samples)
    assert len(samp_ids_to_keep) == len(set(samp_ids_to_keep)), \
           "Duplicate sample IDs in subset"

    return (filter_samples_from_otu_table(otu_table, samp_ids_to_keep, 0, inf),
            filter_mapping_file_from_mapping_f(map_f, samp_ids_to_keep))

def _choose_items_from_bins(sequence, num_items):
    # Adapted from http://stackoverflow.com/a/9873935
    items = []

    for i in range(num_items):
        start = int(ceil(i * float(len(sequence)) / num_items))
        end = int(ceil((i + 1) * float(len(sequence)) / num_items)) - 1
        items.append(sequence[randint(start, end)])

    return items

def generate_simulated_data(sim_data_type, in_dir, out_dir, tests, tree_fp):
    """Simulate gradient or cluster data with simsam.py.

    sim_data_type should be either 'gradient' or 'cluster'.
    """
    create_dir(out_dir)
    otu_table_fp = join(in_dir, tests['study'], 'otu_table.biom')
    map_fp = join(in_dir, tests['study'], 'map.txt')
    map_f = open(map_fp, 'U')
    depth = tests['depth']
    metric = tests['metric']
    num_trials = tests['num_trials']
    category = tests['category']

    # Rarefy the table first since simsam.py's output tables will still have
    # even sampling depth and we don't want to lose simulated samples after the
    # fact.
    even_otu_table_fp = join(out_dir,
            add_filename_suffix(otu_table_fp, '_even%d' % depth))

    if not exists(even_otu_table_fp):
        run_command('single_rarefaction.py -i %s -o %s -d %d;' % (otu_table_fp,
                even_otu_table_fp, depth))

    num_samps = get_num_samples(even_otu_table_fp)

    cmds = []
    for trial_num in range(num_trials):
        trial_num_dir = join(out_dir, '%d' % trial_num)
        create_dir(trial_num_dir)

        for samp_size in tests['sample_sizes']:
            samp_size_dir = join(trial_num_dir, '%d' % samp_size)
            create_dir(samp_size_dir)

            # Lots of duplicate code between these two blocks... need to
            # refactor and test.
            if samp_size <= num_samps:
                simsam_rep_num = 1

                subset_otu_table_fp = join(samp_size_dir,
                                           basename(even_otu_table_fp))
                subset_map_fp = join(samp_size_dir, basename(map_fp))

                if not has_results(samp_size_dir,
                        required_files=[basename(subset_otu_table_fp),
                                        basename(subset_map_fp)]):
                    run_command('choose_data_subset.py -t %s -i %s -m %s '
                                '-c %s -n %d -o %s' % (sim_data_type,
                                                       even_otu_table_fp,
                                                       map_fp, category,
                                                       samp_size,
                                                       samp_size_dir))

                for d in tests['dissim']:
                    dissim_dir = join(samp_size_dir, repr(d))

                    required_files = ['%s_dm.txt' % metric,
                                      '%s_pc.txt' % metric,
                                      'map.txt']
                    if sim_data_type == 'gradient':
                        required_files.append('%s_dm.txt' % category)

                    if not has_results(dissim_dir,
                                       required_files=required_files):
                        simsam_map_fp = join(dissim_dir,
                                add_filename_suffix(subset_map_fp,
                                '_n%d_d%r' % (simsam_rep_num, d)))
                        simsam_otu_table_fp = join(dissim_dir,
                                add_filename_suffix(subset_otu_table_fp,
                                '_n%d_d%r' % (simsam_rep_num, d)))

                        cmd = ('simsam.py -i %s -t %s -o %s -d %r -n %d '
                               '-m %s;' % (subset_otu_table_fp, tree_fp,
                                           dissim_dir, d, simsam_rep_num,
                                           subset_map_fp))

                        if sim_data_type == 'gradient':
                            cmd += ('distance_matrix_from_mapping.py -i %s '
                                    '-c %s -o %s;' % (simsam_map_fp, category,
                                    join(dissim_dir, '%s_dm.txt' % category)))

                        cmd += 'beta_diversity.py -i %s -o %s -m %s -t %s;' % (
                                simsam_otu_table_fp, dissim_dir, metric,
                                tree_fp)
                        cmd += 'mv %s %s;' % (join(dissim_dir, '%s_%s.txt' % (
                                metric,
                                splitext(basename(simsam_otu_table_fp))[0])),
                                '%s_dm.txt' % join(dissim_dir, metric))
                        cmd += 'cp %s %s;' % (simsam_map_fp,
                                             join(dissim_dir, 'map.txt'))
                        cmd += 'principal_coordinates.py -i %s -o %s' % (
                                join(dissim_dir, '%s_dm.txt' % metric),
                                join(dissim_dir, '%s_pc.txt' % metric))

                        cmds.append(cmd)
            else:
                # We need to simulate more samples than we originally have.
                simsam_rep_num = int(ceil(samp_size / num_samps))

                for d in tests['dissim']:
                    dissim_dir = join(samp_size_dir, repr(d))

                    required_files = ['%s_dm.txt' % metric,
                                      '%s_pc.txt' % metric,
                                      'map.txt']
                    if sim_data_type == 'gradient':
                        required_files.append('%s_dm.txt' % category)

                    if not has_results(dissim_dir,
                                       required_files=required_files):
                        simsam_map_fp = join(dissim_dir,
                                add_filename_suffix(map_fp,
                                '_n%d_d%r' % (simsam_rep_num, d)))
                        simsam_otu_table_fp = join(dissim_dir,
                                add_filename_suffix(even_otu_table_fp,
                                '_n%d_d%r' % (simsam_rep_num, d)))

                        cmd = ('simsam.py -i %s -t %s -o %s -d %r -n %d '
                               '-m %s;' % (even_otu_table_fp, tree_fp,
                                           dissim_dir, d, simsam_rep_num,
                                           map_fp))

                        subset_dir = join(dissim_dir, 'subset')
                        cmd += ('choose_data_subset.py -t %s -i %s -m %s '
                                '-c %s -n %d -o %s;' % (sim_data_type,
                                                        simsam_otu_table_fp,
                                                        simsam_map_fp,
                                                        category, samp_size,
                                                        subset_dir))
                        subset_otu_table_fp = join(subset_dir,
                                basename(simsam_otu_table_fp))
                        subset_map_fp = join(subset_dir,
                                basename(simsam_map_fp))

                        if sim_data_type == 'gradient':
                            cmd += ('distance_matrix_from_mapping.py -i %s '
                                    '-c %s -o %s;' % (subset_map_fp, category,
                                                      join(dissim_dir,
                                                      '%s_dm.txt' % category)))

                        cmd += 'beta_diversity.py -i %s -o %s -m %s -t %s;' % (
                                subset_otu_table_fp, dissim_dir, metric,
                                tree_fp)
                        cmd += 'mv %s %s;' % (join(dissim_dir, '%s_%s.txt' % (
                                metric,
                                splitext(basename(subset_otu_table_fp))[0])),
                                '%s_dm.txt' % join(dissim_dir, metric))
                        cmd += 'cp %s %s;' % (subset_map_fp,
                                             join(dissim_dir, 'map.txt'))
                        cmd += 'principal_coordinates.py -i %s -o %s' % (
                                join(dissim_dir, '%s_dm.txt' % metric),
                                join(dissim_dir, '%s_pc.txt' % metric))

                        cmds.append(cmd)

    run_parallel_jobs(cmds, run_command)

def process_simulated_data(in_dir, tests):
    """Run statistical methods over simulated data."""
    metric = tests['metric']
    category = tests['category']
    num_perms = tests['num_perms']
    num_trials = tests['num_trials']

    cmds = []
    for trial_num in range(num_trials):
        trial_num_dir = join(in_dir, '%d' % trial_num)

        for samp_size in tests['sample_sizes']:
            samp_size_dir = join(trial_num_dir, '%d' % samp_size)

            for d in tests['dissim']:
                dissim_dir = join(samp_size_dir, repr(d))

                dm_fp = join(dissim_dir, '%s_dm.txt' % metric)
                map_fp = join(dissim_dir, 'map.txt')
                grad_dm_fp = join(dissim_dir, '%s_dm.txt' % category)

                for method in tests['methods']:
                    method_dir = join(dissim_dir, method)
                    create_dir(method_dir)

                    if not has_results(method_dir):
                        if method == 'mantel' or method == 'mantel_corr':
                            in_dm_fps = ','.join((dm_fp, grad_dm_fp))

                            cmds.append('compare_distance_matrices.py '
                                        '--method %s -n %d -i %s -o %s' % (
                                            method, num_perms, in_dm_fps,
                                            method_dir))
                        elif method == 'partial_mantel' or method == 'best':
                            raise NotImplementedError("%s method is not "
                                                      "currently supported." %
                                                      method)
                        else:
                            cmds.append('compare_categories.py --method %s '
                                        '-i %s -m %s -c %s -o %s -n %d' % (
                                            method, dm_fp, map_fp, category,
                                            method_dir, num_perms))

    run_parallel_jobs(cmds, run_command)

def create_sample_size_plots(sim_data_type, in_dir, tests):
    """Create plots of sample size vs effect size/p-val for each dissim."""
    study = tests['study']
    category = tests['category']

    # +1 to account for PCoA plot.
    num_rows = len(tests['methods']) + 1
    # test stat, p-val, legend
    num_cols = 3

    fig = figure(num=None, figsize=(20, 20), facecolor='w', edgecolor='k')
    fig.suptitle('%s: %s' % (study, category))

    for method_idx, (method, parse_fn) in enumerate(tests['methods'].items()):
        # dissim -> {'sample_sizes': list,
        #            'effect_sizes': list of lists, one for each trial,
        #            'p_vals' -> list of lists, one for each trial}
        plots_data = defaultdict(lambda: defaultdict(list))

        for trial_num in range(tests['num_trials']):
            trial_num_dir = join(in_dir, '%d' % trial_num)

            for samp_size in tests['sample_sizes']:
                samp_size_dir = join(trial_num_dir, '%d' % samp_size)

                for d in tests['dissim']:
                    method_dir = join(samp_size_dir, repr(d), method)

                    effect_size, p_val = parse_fn(open(join(method_dir,
                            '%s_results.txt' % method), 'U'))

                    if samp_size not in plots_data[d]['sample_sizes']:
                        plots_data[d]['sample_sizes'].append(samp_size)
                        plots_data[d]['effect_sizes'].append([])
                        plots_data[d]['p_vals'].append([])

                    samp_size_idx = \
                            plots_data[d]['sample_sizes'].index(samp_size)
                    plots_data[d]['effect_sizes'][samp_size_idx].append(
                            effect_size)
                    plots_data[d]['p_vals'][samp_size_idx].append(p_val)

        # plot_num is 1-based indexing.
        plot_num = method_idx * num_cols + 1
        ax1 = subplot(num_rows, num_cols, plot_num)
        ax2 = subplot(num_rows, num_cols, plot_num + 1)

        color_pool = get_color_pool()

        legend_labels = []
        legend_lines = []
        for d, plot_data in sorted(plots_data.items()):
            avg_effect_sizes = []
            std_effect_sizes = []
            for e in plot_data['effect_sizes']:
                assert len(e) == tests['num_trials']
                avg_effect_sizes.append(mean(e))
                std_effect_sizes.append(std(e))

            avg_p_vals = []
            std_p_vals = []
            for e in plot_data['p_vals']:
                assert len(e) == tests['num_trials']
                avg_p_vals.append(mean(e))
                std_p_vals.append(std(e))

            assert len(plot_data['sample_sizes']) == \
                   len(avg_effect_sizes), "%d != %d" % (
                   len(plot_data['sample_sizes']),
                   len(avg_effect_sizes))

            assert len(plot_data['sample_sizes']) == \
                   len(avg_p_vals), "%d != %d" % (
                   len(plot_data['sample_sizes']),
                   len(avg_p_vals))

            color = color_pool.pop(0)
            label = 'd=%r' % d
            legend_labels.append(label)
            legend_lines.append(Line2D([0, 1], [0, 0], color=color,
                                linewidth=2))

            # Plot test statistics.
            ax1.errorbar(plot_data['sample_sizes'], avg_effect_sizes,
                         yerr=std_effect_sizes, color=color,
                         label=label, fmt='-')

            # Plot p-values.
            _, _, barlinecols = ax2.errorbar(plot_data['sample_sizes'],
                                             avg_p_vals, yerr=std_p_vals,
                                             color=color, label=label,
                                             linestyle='--')
            barlinecols[0].set_linestyles('dashed')

        ax2.set_yscale('log', nonposy='clip')
        x_label = 'Number of samples'
        ax1.set_xlabel(x_label)
        ax2.set_xlabel(x_label)
        ax1.set_ylabel('%s\n\ntest statistic' % method)
        ax2.set_ylabel('p-value')

        if method_idx == 0:
            ax3 = subplot(num_rows, num_cols, plot_num + 2, frame_on=False)
            ax3.get_xaxis().set_visible(False)
            ax3.get_yaxis().set_visible(False)
            ax3.legend(legend_lines, legend_labels, ncol=2, title='Legend',
                       loc='upper left', fancybox=True, shadow=True)

    # Plot PCoA as last row.
    plot_pcoa(sim_data_type, in_dir, tests, num_rows, num_cols)

    fig.tight_layout(pad=5.0, w_pad=2.0, h_pad=2.0)
    fig.savefig(join(in_dir, '%s_%s.pdf' % (tests['study'], category)),
                format='pdf')

def plot_pcoa(sim_data_type, in_dir, tests, num_rows, num_cols):
    plot_num = (num_rows - 1) * num_cols + 1
    trial_num = 0
    samp_size = tests['pcoa_sample_size']
    metric = tests['metric']
    category = tests['category']

    trial_num_dir = join(in_dir, '%d' % trial_num)
    samp_size_dir = join(trial_num_dir, '%d' % samp_size)

    for d_idx, d in enumerate(tests['pcoa_dissim']):
        dissim_dir = join(samp_size_dir, repr(d))
        pc_fp = join(dissim_dir, '%s_pc.txt' % metric)
        map_fp = join(dissim_dir, 'map.txt')

        pc_f = open(pc_fp, 'U')
        map_f = open(map_fp, 'U')
        pc_data = parse_coords(pc_f)
        pc_f.seek(0)

        ax = subplot(num_rows, num_cols, plot_num + d_idx)

        if sim_data_type == 'gradient':
            # Build list of (gradient value, sid) tuples.
            xs, ys, gradient = _collate_gradient_pcoa_plot_data(pc_f, map_f,
                                                                category)
            scatter(xs, ys, s=80, c=gradient, cmap='RdYlBu')
        elif sim_data_type == 'cluster':
            plot_data = _collate_cluster_pcoa_plot_data(pc_f, map_f, category)
            for xs, ys, color in plot_data:
                scatter(xs, ys, color=color)
        else:
            raise ValueError("Unrecognized simulated data type '%s'." %
                             sim_data_type)

        title('d=%r' % d)
        xlabel('PC1 (%1.2f%%)' % pc_data[3][0])
        ylabel('PC2 (%1.2f%%)' % pc_data[3][1])
        xticks([])
        yticks([])

def _collate_gradient_pcoa_plot_data(coords_f, map_f, category):
    pc_data = parse_coords(coords_f)
    coords_d = dict(zip(pc_data[0], pc_data[1]))

    # Build list of (gradient value, sid) tuples.
    map_dict = parse_mapping_file_to_dict(map_f)[0]
    sorted_sids = sorted([(float(md[category]), sid)
                          for sid, md in map_dict.items()])

    xs = [coords_d[sid][0] for _, sid in sorted_sids]
    ys = [coords_d[sid][1] for _, sid in sorted_sids]
    gradient = [cat_val for cat_val, _ in sorted_sids]

    return xs, ys, gradient

def _collate_cluster_pcoa_plot_data(coords_f, map_f, category):
    pc_data = parse_coords(coords_f)
    coords_d = dict(zip(pc_data[0], pc_data[1]))

    map_data = parse_mapping_file(map_f)
    full_map_data = [map_data[1]]
    full_map_data.extend(map_data[0])

    sid_map = group_by_field(full_map_data, category)
    sorted_states = sorted(sid_map.keys())

    color_pool = get_color_pool()
    if len(sorted_states) > len(color_pool):
        raise ValueError("Not enough colors to uniquely color sample "
                         "groups.")

    results = []
    for state, color in zip(sorted_states,
                            color_pool[:len(sorted_states)]):
        sids = sid_map[state]
        xs = [coords_d[sid][0] for sid in sids]
        ys = [coords_d[sid][1] for sid in sids]
        results.append((xs, ys, color))

    return results

def main():
    test = False

    if test:
        in_dir = 'test_datasets'
        out_dir = 'test_simulated_output'
        out_gradient_dir = join(out_dir, 'gradient')
        out_cluster_dir = join(out_dir, 'cluster')
        tree_fp = join('test_datasets', 'overview', 'rep_set.tre')
        gradient_tests = {
            'study': 'overview',
            'depth': 146,
            'metric': 'unweighted_unifrac',
            'num_perms': 999,
            'dissim': [0.0, 0.001, 0.01, 0.1, 1.0, 10.0],
            'pcoa_dissim': [0.001, 0.1, 1.0],
            'sample_sizes': [3, 5, 13],
            'pcoa_sample_size': 13,
            'num_trials': 3,
            'category': 'Gradient',
            'methods': {
                'mantel': parse_mantel_results,
                #'morans_i': parse_morans_i_results
            }
        }

        cluster_tests = {
            'study': 'overview',
            'depth': 146,
            'metric': 'unweighted_unifrac',
            'num_perms': 999,
            'dissim': [0.0, 0.001, 0.01, 0.1, 1.0, 10.0],
            'pcoa_dissim': [0.001, 0.1, 1.0],
            'sample_sizes': [3, 5, 13],
            'pcoa_sample_size': 13,
            'num_trials': 3,
            'category': 'Treatment',
            'methods': {
                'adonis': parse_adonis_results,
                'anosim': parse_anosim_permanova_results
            }
        }
    else:
        in_dir = '../data'
        out_dir = 'sim_data_output'
        out_gradient_dir = join(out_dir, 'gradient')
        out_cluster_dir = join(out_dir, 'cluster')
        tree_fp = join('gg_otus_4feb2011', 'trees', 'gg_97_otus_4feb2011.tre')
        gradient_tests = {
            'study': '88_soils',
            'depth': 400,
            'metric': 'unweighted_unifrac',
            'num_perms': 999,
            # dissim must all be floats!
            'dissim': [0.0, 0.001, 0.01, 0.1, 0.4, 0.7, 1.0, 10.0, 40.0, 70.0,
                       100.0],
            'pcoa_dissim': [0.001, 0.1, 1.0],
            # sample_sizes must all be ints!
            'sample_sizes': [5, 10, 20, 40, 60, 80, 100, 150, 200, 300],
            'pcoa_sample_size': 150,
            'num_trials': 10,
            'category': 'PH',
            'methods': {
                'mantel': parse_mantel_results,
                #'morans_i': parse_morans_i_results
            }
        }

        cluster_tests = {
            'study': 'keyboard',
            'depth': 390,
            'metric': 'unweighted_unifrac',
            'num_perms': 999,
            # dissim must all be floats!
            'dissim': [0.0, 0.001, 0.01, 0.1, 0.4, 0.7, 1.0, 10.0, 40.0, 70.0,
                       100.0],
            'pcoa_dissim': [0.001, 0.1, 1.0],
            # sample_sizes must all be ints!
            'sample_sizes': [5, 10, 20, 40, 60, 80, 100, 150, 200, 300],
            'pcoa_sample_size': 150,
            'num_trials': 10,
            'category': 'HOST_SUBJECT_ID',
            'methods': {
                'adonis': parse_adonis_results,
                'anosim': parse_anosim_permanova_results,
                'mrpp': parse_mrpp_results,
                'permanova': parse_anosim_permanova_results,
                'dbrda': parse_dbrda_results
            }
        }

    #generate_simulated_data('gradient', in_dir, out_gradient_dir,
    #                        gradient_tests, tree_fp)
    #generate_simulated_data('cluster', in_dir, out_cluster_dir, cluster_tests,
    #                        tree_fp)
    #process_simulated_data(out_gradient_dir, gradient_tests)
    #process_simulated_data(out_cluster_dir, cluster_tests)
    create_sample_size_plots('gradient', out_gradient_dir, gradient_tests)
    create_sample_size_plots('cluster', out_cluster_dir, cluster_tests)


if __name__ == "__main__":
    main()
