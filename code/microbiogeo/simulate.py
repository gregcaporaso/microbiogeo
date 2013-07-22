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

from collections import defaultdict
from os import listdir
from os.path import basename, exists, join, splitext
from random import randint, sample

from biom.parse import parse_biom_table

from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.pyplot import figure
from matplotlib.ticker import FormatStrFormatter

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
from microbiogeo.util import (get_color_pool,
                              get_num_samples_in_distance_matrix,
                              get_num_samples_in_map, get_num_samples_in_table,
                              get_panel_label, get_simsam_rep_num, has_results,
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

def create_simulated_data_plots(analysis_type, in_dir, workflow):
    """Create plots of sample size vs effect size/p-val for each dissim.
    
    Plots will be placed directly under in_dir and will be named according to
    the following convention:

    <study>_<category>_<depth>_<metric>.pdf
    """
    for study in workflow:
        study_dir = join(in_dir, study)

        num_trials = workflow[study]['num_sim_data_trials']
        methods = workflow[study]['methods']
        if Best() in methods:
            methods.remove(Best())
        if MantelCorrelogram() in methods:
            methods.remove(MantelCorrelogram())

        num_methods = len(methods)
        num_rows = max(num_methods, len(workflow[study]['pcoa_dissim']) + 1)
        # test stat, p-val, legend/PCoA.
        num_cols = 3

        for depth in workflow[study]['depths']:
            depth_dir = join(study_dir, '%d' % depth[0])
            data_type_dir = join(depth_dir, 'simulated')

            for category in workflow[study]['categories']:
                category_dir = join(data_type_dir, category[0])

                # metric -> Figure
                figs = {}
                for metric in workflow[study]['metrics']:
                    figs[metric[0]] = figure(num=None, figsize=(20, 20),
                                             facecolor='w', edgecolor='k')

                for method_idx, method in enumerate(methods):
                    # metric ->
                    #     dissim -> {
                    #         'sample_sizes': list,
                    #         'effect_sizes': list of lists, one for each size,
                    #         'p_vals' -> list of lists, one for each size
                    #     }
                    plots_data = defaultdict(lambda:
                            defaultdict(lambda: defaultdict(list)))

                    for trial_num in range(num_trials):
                        trial_num_dir = join(category_dir, '%d' % trial_num)

                        for samp_size in workflow[study]['sample_sizes']:
                            samp_size_dir = join(trial_num_dir,
                                                 '%d' % samp_size)

                            for d in workflow[study]['plot_dissim']:
                                dissim_dir = join(samp_size_dir, repr(d))

                                for metric in workflow[study]['metrics']:
                                    metric_dir = join(dissim_dir, metric[0])
                                    method_dir = join(metric_dir, method.Name)

                                    results_fp = join(method_dir,
                                            '%s_results.txt' % method.Name)
                                    effect_size, p_val = method.parse(
                                            open(results_fp, 'U'))

                                    if samp_size not in plots_data[metric[0]][d]['sample_sizes']:
                                        plots_data[metric[0]][d]['sample_sizes'].append(samp_size)
                                        plots_data[metric[0]][d]['effect_sizes'].append([])
                                        plots_data[metric[0]][d]['p_vals'].append([])

                                    samp_size_idx = plots_data[metric[0]][d]['sample_sizes'].index(samp_size)
                                    plots_data[metric[0]][d]['effect_sizes'][samp_size_idx].append(effect_size)
                                    plots_data[metric[0]][d]['p_vals'][samp_size_idx].append(p_val)

                    for metric in workflow[study]['metrics']:
                        fig = figs[metric[0]]
                        metric_plots_data = plots_data[metric[0]]

                        # plot_num is 1-based indexing.
                        plot_num = method_idx * num_cols + 1
                        ax1 = fig.add_subplot(num_rows, num_cols, plot_num)
                        ax2 = fig.add_subplot(num_rows, num_cols, plot_num + 1)

                        color_pool = get_color_pool()

                        min_dissim = min(metric_plots_data.keys())
                        max_dissim = max(metric_plots_data.keys())

                        legend_labels = []
                        legend_lines = []
                        for d, plot_data in sorted(metric_plots_data.items()):
                            avg_effect_sizes, std_effect_sizes, avg_p_vals, std_p_vals = \
                                    _compute_plot_data_statistics(plot_data, num_trials)
                            color = color_pool.pop(0)

                            label = 'd=%r' % d
                            if d == 0.0:
                                label += ' (original data)'
                            #elif d == max_dissim:
                            #    label += ' (neg. control)'

                            legend_labels.append(label)
                            legend_lines.append(Line2D([0, 1], [0, 0],
                                                color=color, linewidth=2))

                            # Make the original data plot a bit thicker than
                            # the rest.
                            if d == 0.0:
                                line_width = 3
                            else:
                                line_width = 1

                            # Plot test statistics.
                            ax1.errorbar(plot_data['sample_sizes'],
                                    avg_effect_sizes, yerr=std_effect_sizes,
                                    color=color, label=label,
                                    linewidth=line_width, fmt='-')

                            # Plot p-values.
                            _, _, barlinecols = ax2.errorbar(
                                    plot_data['sample_sizes'], avg_p_vals,
                                    yerr=std_p_vals, color=color, label=label,
                                    linewidth=line_width, linestyle='--')
                            barlinecols[0].set_linestyles('dashed')

                        ax1.set_xscale('log', nonposx='clip', basex=2)
                        ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
                        ax2.set_xscale('log', nonposx='clip', basex=2)
                        ax2.xaxis.set_major_formatter(FormatStrFormatter('%d'))

                        ax2.set_yscale('log', nonposy='clip')

                        x_label = 'Number of samples'
                        ax1.set_xlabel(x_label)
                        ax2.set_xlabel(x_label)
                        ax1.set_ylabel('%s (%s)' % (method.DisplayName,
                                                    method.StatDisplayName))
                        ax2.set_ylabel('p-value')

                        min_x = min(workflow[study]['sample_sizes'])
                        max_x = max(workflow[study]['sample_sizes'])
                        ax1.set_xlim(min_x - 0.5, max_x)
                        ax2.set_xlim(min_x - 0.5, max_x)

                        for ax_idx, ax in enumerate((ax1, ax2)):
                            panel_idx = method_idx * 2 + ax_idx
                            panel_label = get_panel_label(panel_idx)
                            xmin = ax.get_xlim()[0]
                            ymin, ymax = ax.get_ylim()
                            yrange = ymax - ymin

                            # Not sure why the math isn't working out for the
                            # p-value plots...
                            if ax is ax1:
                                factor = 0.05
                            else:
                                factor = 0.60

                            ax.text(xmin, ymax + (factor * yrange),
                                    '(%s)' % panel_label)

                        if method_idx == 0:
                            ax3 = fig.add_subplot(num_rows, num_cols,
                                                  plot_num + 2, frame_on=False)
                            ax3.get_xaxis().set_visible(False)
                            ax3.get_yaxis().set_visible(False)

                            start_panel_label = get_panel_label(0)
                            end_panel_label = \
                                    get_panel_label(num_methods * 2 - 1)

                            if analysis_type == 'gradient':
                                loc='center'
                            elif analysis_type == 'cluster':
                                loc='center left'

                            assert len(legend_lines) == len(workflow[study]['plot_dissim'])
                            assert len(legend_labels) == len(workflow[study]['plot_dissim'])
                            legend_title = ('           Legend (Panels %s-%s)\nd = '
                                    '"noise" introduced to samples' % (
                                        start_panel_label, end_panel_label))
                            ax3.legend(legend_lines, legend_labels, ncol=1,
                                    title=legend_title, loc=loc, fancybox=True,
                                    shadow=True)

                for metric in workflow[study]['metrics']:
                    fig = figs[metric[0]]

                    # Plot PCoA in last column of figure.
                    plot_pcoa(analysis_type, fig, category_dir, workflow[study],
                            category, metric, num_rows, num_cols, num_methods)

                    fig.tight_layout(pad=5.0, w_pad=2.0, h_pad=2.0)
                    fig.savefig(join(in_dir, '%s_%s_%d_%s.pdf' % (study,
                            category[0], depth[0], metric[0])), format='pdf')
                    fig.savefig(join(in_dir, '%s_%s_%d_%s.png' % (study,
                            category[0], depth[0], metric[0])), format='png',
                            dpi=100)

def _compute_plot_data_statistics(plot_data, num_trials):
    avg_effect_sizes = []
    std_effect_sizes = []
    for e in plot_data['effect_sizes']:
        if len(e) != num_trials:
            raise ValueError("Data length doesn't match the number of trials.")

        avg_effect_sizes.append(mean(e))
        std_effect_sizes.append(std(e))

    # Need to compute asymmetric error bars for p-values to avoid
    # negative error bars on log scale.
    avg_p_vals = []
    std_p_vals = [[], []]
    for e in plot_data['p_vals']:
        if len(e) != num_trials:
            raise ValueError("Data length doesn't match the number of trials.")

        avg_p_val = mean(e)
        std_p_val = std(e)
        lower_std_p_val = std_p_val
        upper_std_p_val = std_p_val

        # Cut off lower bound at 1e-5. This is important for some methods such
        # as Moran's I where the p-value can get *very* close to zero since it
        # isn't based on permutations.
        if avg_p_val < 1e-5:
            avg_p_val = 1e-5
            lower_std_p_val = 0.0
            upper_std_p_val = 0.0
        elif avg_p_val - std_p_val < 1e-5:
            lower_std_p_val = avg_p_val - 1e-5

        avg_p_vals.append(avg_p_val)
        std_p_vals[0].append(lower_std_p_val)
        std_p_vals[1].append(upper_std_p_val)

    if len(plot_data['sample_sizes']) != len(avg_effect_sizes):
        raise ValueError("%d != %d" % (len(plot_data['sample_sizes']),
                                       len(avg_effect_sizes)))

    if len(plot_data['sample_sizes']) != len(avg_p_vals):
        raise ValueError("%d != %d" % (len(plot_data['sample_sizes']),
                                       len(avg_p_vals)))

    return avg_effect_sizes, std_effect_sizes, avg_p_vals, std_p_vals

def plot_pcoa(analysis_type, fig, in_dir, workflow, category, metric, num_rows,
              num_cols, num_methods):
    trial_num = 0
    samp_size = workflow['pcoa_sample_size']

    trial_num_dir = join(in_dir, '%d' % trial_num)
    samp_size_dir = join(trial_num_dir, '%d' % samp_size)

    min_dissim = min(workflow['pcoa_dissim'])
    max_dissim = max(workflow['pcoa_dissim'])

    legend_symbols = []
    legend_labels = []
    for d_idx, d in enumerate(workflow['pcoa_dissim']):
        dissim_dir = join(samp_size_dir, repr(d))
        metric_dir = join(dissim_dir, metric[0])

        pc_fp = join(metric_dir, 'pc.txt')
        map_fp = join(metric_dir, 'map.txt')

        pc_f = open(pc_fp, 'U')
        map_f = open(map_fp, 'U')
        pc_data = parse_coords(pc_f)
        pc_f.seek(0)
        assert len(pc_data[0]) == samp_size

        # Skip the first row (the legend is already at that cell).
        plot_num = (d_idx + 2) * num_cols
        ax = fig.add_subplot(num_rows, num_cols, plot_num)

        if analysis_type == 'gradient':
            # Build list of (gradient value, sid) tuples.
            xs, ys, gradient = _collate_gradient_pcoa_plot_data(pc_f, map_f,
                                                                category[0])
            scatter_colorbar_data = ax.scatter(xs, ys, s=80, c=gradient,
                                               cmap='RdYlBu')
            # We have to use gridspec to get this to work with tight_layout.
            cb = fig.colorbar(scatter_colorbar_data, use_gridspec=True)
            cb.set_label(category[1])
        elif analysis_type == 'cluster':
            plot_data = _collate_cluster_pcoa_plot_data(pc_f, map_f,
                                                        category[0])
            for xs, ys, color, state in plot_data:
                ax.scatter(xs, ys, color=color, label=state)

                if d_idx == 0:
                    legend_symbols.append(Line2D(range(1), range(1),
                                          color='white', marker='o',
                                          markeredgecolor=color,
                                          markerfacecolor=color))
                    legend_labels.append(category[2].get(state, state))
        else:
            raise ValueError("Unrecognized simulated data type '%s'." %
                             analysis_type)

        plot_title = 'd=%r' % d
        if d == 0.0:
            plot_title += ' (original data)'
        #elif d == max_dissim:
        #    plot_title += ' (neg. control)'
        ax.set_title(plot_title)

        ax.set_xlabel('PC1 (%1.2f%%)' % pc_data[3][0])
        ax.set_ylabel('PC2 (%1.2f%%)' % pc_data[3][1])
        ax.set_xticks([])
        ax.set_yticks([])

        panel_idx = num_methods * 2 + d_idx
        panel_label = get_panel_label(panel_idx)
        xmin = ax.get_xlim()[0]
        ymin, ymax = ax.get_ylim()
        yrange = ymax - ymin
        ax.text(xmin, ymax + (0.04 * yrange), '(%s)' % panel_label)

    if analysis_type == 'cluster':
        # Plot our new legend and add the existing one back.
        legend_ax = fig.add_subplot(num_rows, num_cols, 3, frame_on=False)
        existing_legend = legend_ax.get_legend()
        existing_legend.set_bbox_to_anchor((-0.05, 0.5))

        start_panel_label = get_panel_label(num_methods * 2)
        end_panel_label = get_panel_label(num_methods * 2 +
                                          len(workflow['pcoa_dissim']) - 1)

        assert len(legend_symbols) == len(legend_labels)
        legend_ax.legend(legend_symbols, legend_labels, ncol=1,
                   title='Legend (Panels %s-%s)' % (start_panel_label,
                                                    end_panel_label),
                   loc='center right', fancybox=True, shadow=True, numpoints=1,
                   bbox_to_anchor=(1.05, 0.5))

        legend_ax.add_artist(existing_legend)

    # Draw box around PCoA plots. Do the math in figure coordinates.
    top_ax = fig.add_subplot(num_rows, num_cols, 6)
    rec = Rectangle((1 - (1 / num_cols) + 0.005, 0),
                    (1 / num_cols) - 0.005,
                    1 - (1 / num_rows) - 0.005,
                    fill=False, lw=2, clip_on=False,
                    transform=top_ax.figure.transFigure)
    top_ax.add_patch(rec)

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
        results.append((xs, ys, color, state))

    return results
