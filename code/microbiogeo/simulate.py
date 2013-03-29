#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Module with functionality for simulating data."""

from biom.parse import parse_biom_table
from collections import defaultdict
from numpy import ceil, inf
from os import listdir
from os.path import basename, join, splitext
from matplotlib.pyplot import figure, legend, title, xlim
from qiime.colors import data_colors, data_color_order
from qiime.filter import (filter_mapping_file_from_mapping_f,
                          filter_samples_from_otu_table)
from qiime.make_distance_histograms import matplotlib_rgb_color
from qiime.parse import parse_mapping_file_to_dict
from qiime.util import add_filename_suffix, create_dir, MetadataMap
from random import randint, sample

from microbiogeo.parse import parse_mantel_results, parse_morans_i_results
from microbiogeo.util import has_results, run_command, run_parallel_jobs

def choose_cluster_subsets(otu_table_f, map_f, category,
                           num_samples_per_group):
    otu_table = parse_biom_table(otu_table_f)
    metadata_map = MetadataMap.parseMetadataMap(map_f)

    category_map = defaultdict(list)
    for samp_id in metadata_map.SampleIds:
        # Mapping files can have more samples than OTU tables.
        if samp_id in otu_table.SampleIds:
            category_val = metadata_map.getCategoryValue(samp_id, category)
            category_map[category_val].append(samp_id)

    samp_ids_to_keep = []
    for category_val, samp_ids in category_map.items():
        samp_ids_to_keep.extend(
                sample(samp_ids, min(num_samples_per_group, len(samp_ids))))

    return (filter_samples_from_otu_table(otu_table, samp_ids_to_keep, 0, inf),
            filter_mapping_file_from_mapping_f(map_f, samp_ids_to_keep),
            len(samp_ids_to_keep))

def choose_gradient_subset(otu_table_f, map_f, category, num_total_samples):
    otu_table = parse_biom_table(otu_table_f)
    mdm, _ = parse_mapping_file_to_dict(map_f)

    try:
        map_f.seek(0)
    except AttributeError:
        pass

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

def generate_gradient_simulated_data(in_dir, out_dir, tests, tree_fp):
    """Simulate gradient data with simsam.py."""
    create_dir(out_dir)
    otu_table_fp = join(in_dir, tests['study'], 'otu_table.biom')
    map_fp = join(in_dir, tests['study'], 'map.txt')
    map_f = open(map_fp, 'U')
    depth = tests['depth']
    metric = tests['metric']
    category = tests['category']

    # Rarefy the table first since simsam.py's output tables will still have
    # even sampling depth and we don't want to lose simulated samples after the
    # fact.
    even_otu_table_fp = join(out_dir,
            add_filename_suffix(otu_table_fp, '_even%d' % depth))
    run_command('single_rarefaction.py -i %s -o %s -d %d;' % (otu_table_fp,
            even_otu_table_fp, depth))
    
    # Figure out how many samples we have in the rarefied table.
    even_otu_table_f = open(even_otu_table_fp, 'U')
    even_otu_table = parse_biom_table(even_otu_table_f)
    even_otu_table_f.seek(0)
    num_samps = len(even_otu_table.SampleIds)

    cmds = []
    for samp_size in tests['sample_sizes']:
        samp_size_dir = join(out_dir, '%d' % samp_size)
        create_dir(samp_size_dir)

        # Lots of duplicate code between these two blocks... need to refactor
        # and test.
        if samp_size <= num_samps:
            simsam_rep_num = 1

            run_command('choose_gradient_subset.py -i %s -m %s -c %s -n %d '
                        '-o %s' % (even_otu_table_fp, map_fp, category,
                                   samp_size, samp_size_dir))
            subset_otu_table_fp = join(samp_size_dir,
                                       basename(even_otu_table_fp))
            subset_map_fp = join(samp_size_dir, basename(map_fp))

            for d in tests['dissim']:
                dissim_dir = join(samp_size_dir, repr(d))
                simsam_map_fp = join(dissim_dir,
                        add_filename_suffix(subset_map_fp,
                                            '_n%d_d%r' % (simsam_rep_num, d)))
                simsam_otu_table_fp = join(dissim_dir,
                        add_filename_suffix(subset_otu_table_fp,
                                            '_n%d_d%r' % (simsam_rep_num, d)))

                cmd = 'simsam.py -i %s -t %s -o %s -d %r -n %d -m %s;' % (
                        subset_otu_table_fp, tree_fp, dissim_dir, d,
                        simsam_rep_num, subset_map_fp)
                cmd += 'distance_matrix_from_mapping.py -i %s -c %s -o %s;' % (
                        simsam_map_fp, category,
                        join(dissim_dir, '%s_dm.txt' % category))
                cmd += 'beta_diversity.py -i %s -o %s -m %s -t %s;' % (
                        simsam_otu_table_fp, dissim_dir, metric, tree_fp)
                cmd += 'mv %s %s;' % (join(dissim_dir, '%s_%s.txt' % (
                        metric, splitext(basename(simsam_otu_table_fp))[0])),
                        '%s_dm.txt' % join(dissim_dir, metric))
                cmd += 'cp %s %s' % (simsam_map_fp,
                                     join(dissim_dir, 'map.txt'))
                cmds.append(cmd)
        else:
            # We need to simulate more samples than we originally have.
            simsam_rep_num = int(ceil(samp_size / num_samps))

            for d in tests['dissim']:
                dissim_dir = join(samp_size_dir, repr(d))
                simsam_map_fp = join(dissim_dir, add_filename_suffix(map_fp,
                        '_n%d_d%r' % (simsam_rep_num, d)))
                simsam_otu_table_fp = join(dissim_dir,
                        add_filename_suffix(even_otu_table_fp,
                                            '_n%d_d%r' % (simsam_rep_num, d)))

                cmd = 'simsam.py -i %s -t %s -o %s -d %r -n %d -m %s;' % (
                        even_otu_table_fp, tree_fp, dissim_dir, d,
                        simsam_rep_num, map_fp)

                subset_dir = join(dissim_dir, 'subset')
                cmd += ('choose_gradient_subset.py -i %s -m %s -c %s -n %d '
                        '-o %s;' % (simsam_otu_table_fp, simsam_map_fp,
                                    category, samp_size, subset_dir))
                subset_otu_table_fp = join(subset_dir,
                                           basename(simsam_otu_table_fp))
                subset_map_fp = join(subset_dir, basename(simsam_map_fp))

                cmd += 'distance_matrix_from_mapping.py -i %s -c %s -o %s;' % (
                        subset_map_fp, category,
                        join(dissim_dir, '%s_dm.txt' % category))
                cmd += 'beta_diversity.py -i %s -o %s -m %s -t %s;' % (
                        subset_otu_table_fp, dissim_dir, metric, tree_fp)
                cmd += 'mv %s %s;' % (join(dissim_dir, '%s_%s.txt' % (
                        metric, splitext(basename(subset_otu_table_fp))[0])),
                        '%s_dm.txt' % join(dissim_dir, metric))
                cmd += 'cp %s %s' % (subset_map_fp,
                                     join(dissim_dir, 'map.txt'))
                cmds.append(cmd)

    run_parallel_jobs(cmds, run_command)

def process_gradient_simulated_data(in_dir, tests):
    """Run statistical methods over gradient simulated data."""
    metric = tests['metric']
    category = tests['category']
    num_perms = tests['num_perms']

    cmds = []
    for samp_size in tests['sample_sizes']:
        samp_size_dir = join(in_dir, '%d' % samp_size)

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

                        cmds.append('compare_distance_matrices.py --method %s '
                                    '-n %d -i %s -o %s' % (method, num_perms,
                                                           in_dm_fps,
                                                           method_dir))
                    elif method == 'morans_i':
                        cmds.append('compare_categories.py --method %s -i %s '
                                    '-m %s -c %s -o %s' % (method, dm_fp,
                                                           map_fp, category,
                                                           method_dir))

    run_parallel_jobs(cmds, run_command)

def create_sample_size_plots(in_dir, tests):
    """Create plots of sample size vs effect size/p-val for each dissim."""
    category = tests['category']

    for method, parse_fn in tests['methods'].items():
        # dissim -> {'sample_sizes': list,
        #            'effect_sizes': list,
        #            'p_vals' -> list}
        plots_data = defaultdict(lambda: defaultdict(list))

        for samp_size in tests['sample_sizes']:
            samp_size_dir = join(in_dir, '%d' % samp_size)

            for d in tests['dissim']:
                method_dir = join(samp_size_dir, repr(d), method)

                effect_size, p_val = parse_fn(
                        open(join(method_dir, '%s_results.txt' % method), 'U'))
                plots_data[d]['sample_sizes'].append(samp_size)
                plots_data[d]['effect_sizes'].append(effect_size)
                plots_data[d]['p_vals'].append(p_val)

        # Twin y-axis code is based on
        # http://matplotlib.org/examples/api/two_scales.html
        fig = figure()
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twinx()

        color_pool = [matplotlib_rgb_color(data_colors[color].toRGB())
                      for color in data_color_order]

        for d, plot_data in sorted(plots_data.items()):
            color = color_pool.pop(0)

            # Plot test statistics on left axis.
            ax1.errorbar(plot_data['sample_sizes'], plot_data['effect_sizes'],
                         color=color, label='d=%r' % d, fmt='-')

            # Plot p-values on the right axis.
            ax2.errorbar(plot_data['sample_sizes'], plot_data['p_vals'],
                         color=color, label='d=%r' % d, linestyle='--')

        #xlim(0, max(plot_data['sample_sizes']))
        #ax2.set_ylim(0.0, 1.0)
        title('%s: %s: %s' % (tests['study'], method, category))
        #lines, labels = ax1.get_legend_handles_labels()
        #lines2, labels2 = ax2.get_legend_handles_labels()
        #ax2.legend(lines + lines2, labels + labels2)
        ax1.set_xlabel('Number of samples')
        ax1.set_ylabel('test statistic')
        ax2.set_ylabel('p-value')
        legend()
        fig.savefig(join(in_dir, '%s_%s_%s.pdf' % (tests['study'], method,
                                                   category)), format='pdf')

def main():
    test = False

    if test:
        in_dir = 'test_datasets'
        out_dir = 'test_simulated_output'
        tree_fp = join('test_datasets', 'overview', 'rep_set.tre')
        gradient_tests = {
            'study': 'overview',
            'depth': 146,
            'metric': 'unweighted_unifrac',
            'num_perms': 999,
            'dissim': [0.001, 0.01, 0.1],
            'sample_sizes': [3, 5, 13, 100],
            'category': 'Gradient',
            'methods': {
                'mantel': parse_mantel_results,
                'morans_i': parse_morans_i_results
            }
        }
    else:
        in_dir = '../datasets'
        out_dir = 'sim_data_output'
        tree_fp = join('gg_otus_4feb2011', 'trees', 'gg_97_otus_4feb2011.tre')
        gradient_tests = {
            'study': '88_soils',
            'depth': 400,
            'metric': 'unweighted_unifrac',
            'num_perms': 999,
            'dissim': [0.001, 0.01, 0.1, 0.5, 1, 10],
            'sample_sizes': [5, 10, 20, 40, 60, 80, 100, 150, 200, 300],
            'category': 'PH',
            'methods': {
                'mantel': parse_mantel_results,
                'morans_i': parse_morans_i_results
            }
        }

    generate_gradient_simulated_data(in_dir, out_dir, gradient_tests, tree_fp)
    process_gradient_simulated_data(out_dir, gradient_tests)
    create_sample_size_plots(out_dir, gradient_tests)


if __name__ == "__main__":
    main()
