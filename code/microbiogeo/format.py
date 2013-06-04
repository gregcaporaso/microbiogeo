#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Module to format data for presentation."""

from collections import defaultdict
from csv import writer
from os.path import join

from numpy import invert, mean, ones, std, tri
import numpy.ma

from matplotlib import cm
from matplotlib.pyplot import (colorbar, figure, imshow, legend, matshow,
                               savefig, subplot, tight_layout, title, xlim,
                               xticks, yticks)

from cogent.maths.stats.test import pearson, spearman

from qiime.util import create_dir

from microbiogeo.util import is_empty

# Not unit-tested.
def create_results_summary_tables(results, out_dir):
    """Creates tables to summarize the results of the statistical methods.

    These tables will be in TSV format so that they can be easily imported into
    Excel for viewing and cleanup for publication.

    A table will be created for each sampling depth / metric combination, for
    each method type (grouping or gradient).
    """
    for depth_desc, depth_res in results.items():
        for metric, metric_res in depth_res.items():
            for method_type, method_type_res in metric_res.items():
                table_rows = format_method_comparison_table(method_type_res)

                table_name = ('%s_analysis_method_comparison_table_%s_%s.txt' %
                              (method_type, depth_desc, metric))
                with open(join(out_dir, table_name), 'wb') as out_f:
                    # We use \r so that we can force linebreaks within cells
                    # when imported into Excel. Not sure if this will work with
                    # other spreadsheet programs such as Open Office.
                    tsv_writer = writer(out_f, delimiter='\t',
                                        lineterminator='\r')
                    tsv_writer.writerows(table_rows)

def format_method_comparison_table(per_method_results):
    constructed_header = None
    rows = []

    for method, method_res in sorted(per_method_results.items()):
        header = ['Method']
        row = [method]

        for study, study_res in sorted(method_res.items()):
            for category, category_res in sorted(study_res.items()):
                header_column_title = '%s\r%s' % (study, category)
                header.append(header_column_title)
                header.append(header_column_title + ' (shuffled)')
                header.append(header_column_title + ' (subsampled)')

                if len(category_res) == 0:
                    row.extend(['N/A'] * 3)
                else:
                    if category_res['full'].isEmpty():
                        row.append('N/A')
                    else:
                        row.append(str(category_res['full']))

                    if category_res['shuffled'].isEmpty():
                        row.append('N/A')
                    else:
                        row.append(str(category_res['shuffled']))

                    cells = []
                    for res in category_res['subsampled']:
                        if res.isEmpty():
                            cells.append('N/A')
                        else:
                            cells.append(str(res))
                    row.append('\r'.join(cells))

        if constructed_header is None:
            rows.append(header)
            constructed_header = header
        elif constructed_header != header:
            raise ValueError("The studies and/or categories did not match up "
                             "exactly between one or more of the methods.")

        rows.append(row)

    return rows

# Not unit-tested.
def create_method_comparison_heatmaps(results, methods, out_dir):
    """Generates heatmaps showing the correlation between each pair of methods.

    Generates two heatmaps (one for Pearson correlation, one for Spearman
    correlation). Uses all available results (e.g. all even sampling depths,
    metrics, and datasets) that match between each pair of methods as input to
    the correlation coefficient methods.

    A heatmap will be written to out_dir for each type of method (grouping or
    gradient).
    """
    for method_type, data in \
            format_method_comparison_heatmaps(results, methods).items():
        for correlation_method, heatmap_data in data.items():
            # Generate the heatmap. Code based on
            # http://matplotlib.org/users/tight_layout_guide.html and
            # http://psaffrey.wordpress.com/2010/07/05/chromosome-interactions-
            #   heatmaps-and-matplotlib/
            fig = figure()
            ax = subplot(111)
            cmap = cm.get_cmap()
            cmap.set_bad('w') # default value is 'k'
            im = ax.imshow(heatmap_data, cmap=cmap, interpolation='nearest')
            method_labels = [method.DisplayName
                             for method in methods[method_type]]

            colorbar(im, use_gridspec=True)
            xticks(range(len(method_labels)), method_labels, rotation=90)
            yticks(range(len(method_labels)), method_labels)

            for loc, spine in ax.spines.items():
                if loc in ['right','top']:
                    spine.set_color('none') # don't draw spine

            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            ax.grid(True, which='minor')

            tight_layout()
            savefig(join(out_dir, '%s_analysis_heatmap_%s.pdf' % (method_type,
                    correlation_method)), format='pdf')
            savefig(join(out_dir, '%s_analysis_heatmap_%s.png' % (method_type,
                    correlation_method)), format='png', dpi=1000)

def format_method_comparison_heatmaps(results, methods):
    shared_studies = {}
    shared_categories = {}

    for depth_desc, depth_res in results.items():
        for metric, metric_res in depth_res.items():
            for method_type, method_type_res in metric_res.items():
                if method_type not in shared_categories:
                    shared_categories[method_type] = {}

                for method, method_res in method_type_res.items():
                    matched_method = False
                    for m in methods[method_type]:
                        if method == m.Name:
                            matched_method = True
                            break

                    if matched_method:
                        studies = sorted(method_res.keys())

                        if method_type not in shared_studies:
                            shared_studies[method_type] = studies
                        elif studies != shared_studies[method_type]:
                            raise ValueError("Not all methods to include in "
                                             "the heatmap have results for "
                                             "the same studies.")

                        for study, study_res in sorted(method_res.items()):
                            categories = [cat for cat, cat_res in \
                                          sorted(study_res.items()) if not
                                          is_empty(cat_res)]

                            if study not in shared_categories[method_type]:
                                shared_categories[method_type][study] = \
                                        set(categories)
                            else:
                                shared_categories[method_type][study] &= \
                                        set(categories)

    # Gather all test statistics for each method (in the same order for each
    # method!).
    method_data = defaultdict(lambda: defaultdict(list))
    for depth_desc, depth_res in results.items():
        for metric, metric_res in depth_res.items():
            for method_type, method_type_res in metric_res.items():
                for method, method_res in method_type_res.items():
                    matched_method = False
                    for m in methods[method_type]:
                        if method == m.Name:
                            matched_method = True
                            break

                    if not matched_method:
                        continue

                    for study, study_res in sorted(method_res.items()):
                        for category, category_res in \
                                sorted(study_res.items()):
                            if category in \
                                    shared_categories[method_type][study]:
                                method_data[method_type][method].append(
                                        category_res['full'].effect_size)
                                method_data[method_type][method].append(
                                        category_res['shuffled'].effect_size)

                                for res in category_res['subsampled']:
                                    method_data[method_type][method].append(
                                            res.effect_size)

    # Make sure our data looks sane. We should have the same number of
    # observations for each method.
    for method_type, results in method_data.items():
        data_length = None

        for method, data in results.items():
            if data_length is None:
                data_length = len(data)
            elif len(data) != data_length:
                raise ValueError("The number of test statistics is not the "
                                 "same between all methods, so we can't "
                                 "compare them.")

    # Compute the correlation coefficient between each pair of methods and put
    # the output in an array. This array can then be used to generate a
    # text-based table or heatmap.
    heatmaps = {}
    for method_type in methods:
        heatmaps[method_type] = {}

        for correlation_name, correlation_fn in \
                ('pearson', pearson), ('spearman', spearman):
            num_methods = len(methods[method_type])
            heatmap_data = ones((num_methods, num_methods))

            # I know this is inefficient, but it really doesn't matter for what
            # we're doing here.
            for method1_idx, method1 in enumerate(methods[method_type]):
                for method2_idx, method2 in enumerate(methods[method_type]):
                    corr_coeff = correlation_fn(
                            method_data[method_type][method1.Name],
                            method_data[method_type][method2.Name])
                    heatmap_data[method1_idx][method2_idx] = corr_coeff

            # Mask out the upper triangle. Taken from
            # http://stackoverflow.com/a/2332520
            mask = invert(tri(heatmap_data.shape[0], k=0, dtype=bool))
            heatmap_data = numpy.ma.array(heatmap_data, mask=mask)

            heatmaps[method_type][correlation_name] = heatmap_data

    return heatmaps

def create_sample_size_plots(in_dir, out_dir, sample_size_tests):
    create_dir(out_dir)

    for method_type, method_type_tests in sample_size_tests.items():
        study = method_type_tests['study']
        metric = method_type_tests['metric']
        categories = method_type_tests['categories']
        subset_sizes = method_type_tests['subset_sizes']
        num_subsets = method_type_tests['num_subsets']
        methods = method_type_tests['methods']

        out_method_type_dir = join(out_dir, method_type)
        out_study_dir = join(out_method_type_dir, study)
        create_dir(out_method_type_dir)
        create_dir(out_study_dir)

        for method in methods:
            # Twin y-axis code is based on
            # http://matplotlib.org/examples/api/two_scales.html
            fig = figure()
            ax1 = fig.add_subplot(111)
            ax2 = ax1.twinx()

            for category, plot_options in categories.items():
                avg_test_stats = []
                std_test_stats = []
                avg_p_vals = []
                std_p_vals = []

                for subset_size in subset_sizes:
                    test_stats = []
                    p_vals = []

                    for subset_num in range(1, num_subsets + 1):
                        if method_type == 'grouping':
                            results_dir = join(in_dir, method_type, study,
                                    '%s_dm_%s_gs%d_%d_%s' % (metric,
                                                             category,
                                                             subset_size,
                                                             subset_num,
                                                             method.Name))
                        elif method_type == 'gradient':
                            results_dir = join(in_dir, method_type, study,
                                    '%s_dm_%s_n%d_%d_%s' % (metric,
                                                            category,
                                                            subset_size,
                                                            subset_num,
                                                            method.Name))
                        else:
                            raise ValueError("Unknown method type '%s'." %
                                             method_type)

                        test_stat, p_val = method.parse(open(join(results_dir,
                                '%s_results.txt' % method.Name), 'U'))
                        test_stats.append(test_stat)
                        p_vals.append(p_val)

                    avg_test_stats.append(mean(test_stats))
                    std_test_stats.append(std(test_stats))
                    avg_p_vals.append(mean(p_vals))
                    std_p_vals.append(std(p_vals))

                # Plot test statistics on left axis.
                ax1.errorbar(subset_sizes, avg_test_stats, yerr=std_test_stats,
                             color=plot_options[0], label=plot_options[1],
                             fmt='-')

                # Plot p-values on the right axis.
                ax2.errorbar(subset_sizes, avg_p_vals, yerr=std_p_vals,
                             color=plot_options[0], label=plot_options[1],
                             fmt='-', linestyle='--')

            xlim(0, 85)
            #ax2.set_ylim(0.0, 1.0)
            title('%s: %s' % (study, method.DisplayName))
            #lines, labels = ax1.get_legend_handles_labels()
            #lines2, labels2 = ax2.get_legend_handles_labels()
            #ax2.legend(lines + lines2, labels + labels2)

            if method_type == 'grouping':
                x_label = 'Samples per group'
            elif method_type == 'gradient':
                x_label = 'Number of samples'
            else:
                raise ValueError("Unknown method type '%s'." % method_type)

            ax1.set_xlabel(x_label)
            ax1.set_ylabel('Average test statistic with standard deviation')
            ax2.set_ylabel('Average p-value with standard deviation')
            legend()
            fig.savefig(join(out_study_dir, '%s_analysis_plot_%s_%s.pdf' % (
                    method_type, study, method.Name)), format='pdf')
