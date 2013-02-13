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

from numpy import ones

from matplotlib.pyplot import (colorbar, figure, imshow, matshow, savefig,
                               subplot, tight_layout, xticks, yticks)

from cogent.maths.stats.test import pearson, spearman

# Not unit-tested.
def create_results_summary_tables(results, out_dir, filename_prefix):
    """Creates tables to summarize the results of the statistical methods.

    These tables will be in TSV format so that they can be easily imported into
    Excel for viewing and cleanup for publication.

    A table will be created for each sampling depth / metric combination and
    written to out_dir with the filename convention
    <filename_prefix>_<depth>_<metric>.txt.
    """
    for depth_desc, depth_res in results.items():
        for metric, metric_res in depth_res.items():
            table_rows = format_method_comparison_table(metric_res)

            with open(join(out_dir, '%s_%s_%s.txt' % (filename_prefix,
                                                      depth_desc,
                                                      metric)), 'wb') as out_f:
                # We use \r so that we can force linebreaks within cells when
                # imported into Excel. Not sure if this will work with other
                # spreadsheet programs such as Open Office.
                tsv_writer = writer(out_f, delimiter='\t', lineterminator='\r')
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

                if len(category_res) == 3:
                    row.append(str(category_res['full']))
                    row.append(str(category_res['shuffled']))
                    row.append('\r'.join(map(str, category_res['subsampled'])))
                else:
                    row.append(['N/A'] * 3)

        if constructed_header is None:
            rows.append(header)
            constructed_header = header
        elif constructed_header != header:
            raise ValueError("The studies and/or categories did not match up "
                             "exactly between one or more of the methods.")

        rows.append(row)

    return rows

def format_p_value_as_asterisk(p_value):
    if not isinstance(p_value, float):
        raise TypeError("p-value must be a float.")
    if p_value < 0 or p_value > 1:
        raise ValueError("p-value must be a float between 0 and 1, inclusive.")

    result = 'x'

    if p_value <= 0.1:
        result = '*'
    if p_value <= 0.05:
        result += '*'
    if p_value <= 0.01:
        result += '*'
    if p_value <= 0.001:
        result += '*'

    return result

# Not unit-tested.
def create_method_comparison_heatmaps(results, methods, method_labels, out_dir,
                                      filename_prefix):
    """Generates heatmaps showing the correlation between each pair of methods.

    Generates two heatmaps (one for Pearson correlation, one for Spearman
    correlation). Uses all available results (e.g. all even sampling depths,
    metrics, and datasets) that match between each pair of methods as input to
    the correlation coefficient methods.

    Each heatmap will be written to out_dir with the filename convention
    <filename_prefix>_<correlation_method>.pdf.
    """
    for correlation_method, heatmap_data in \
            format_method_comparison_heatmaps(results, methods):
        # Generate the heatmap. Code based on
        # http://matplotlib.org/users/tight_layout_guide.html and
        # http://psaffrey.wordpress.com/2010/07/05/chromosome-interactions-
        #   heatmaps-and-matplotlib/
        fig = figure()
        ax = subplot(111)
        im = ax.matshow(heatmap_data, vmin=0, vmax=1)

        colorbar(im, use_gridspec=True)
        xticks(range(len(method_labels)), method_labels, rotation=90)
        yticks(range(len(method_labels)), method_labels)

        tight_layout()
        savefig('%s_%s.pdf' % (filename_prefix, correlation_method),
                format='pdf')

def format_method_comparison_heatmaps(results, methods):
    shared_studies = None
    shared_categories = {}

    for depth_desc, depth_res in results.items():
        for metric, metric_res in depth_res.items():
            for method, method_res in metric_res.items():
                if method in methods:
                    studies = sorted(method_res.keys())

                    if shared_studies is None:
                        shared_studies = studies
                    elif studies != shared_studies:
                        raise ValueError("Not all methods to include in the "
                                         "heatmap have results for the same "
                                         "studies.")

                    for study, study_res in sorted(method_res.items()):
                        categories = [category for category, category_res in
                                      sorted(study_res.items())
                                      if len(category_res) > 0]
                        if study not in shared_categories:
                            shared_categories[study] = set(categories)
                        else:
                            shared_categories[study] &= set(categories)

    # Gather all test statistics for each method (in the same order for each
    # method!).
    method_data = defaultdict(list)
    for depth_desc, depth_res in results.items():
        for metric, metric_res in depth_res.items():
            for method, method_res in metric_res.items():
                if method in methods:
                    for study, study_res in sorted(method_res.items()):
                        for category, category_res in sorted(
                                study_res.items()):
                            if category in shared_categories[study]:
                                method_data[method].append(
                                        category_res['full'].effect_size)
                                method_data[method].append(
                                        category_res['shuffled'].effect_size)
                                
                                for res in category_res['subsampled']:
                                    method_data[method].append(res.effect_size)

    # Make sure our data looks sane. We should have the same number of
    # observations for each method.
    data_length = None
    for method, data in method_data.items():
        if data_length is None:
            data_length = len(data)
        elif len(data) != data_length:
            raise ValueError("The number of test statistics is not the same "
                             "between all methods, so we can't compare them.")

    # Compute the correlation coefficient between each pair of methods and put
    # the output in an array. This array can then be used to generate a
    # text-based table or heatmap.
    heatmaps = {}
    for correlation_name, correlation_fn in \
            ('pearson', pearson), ('spearman', spearman):
        num_methods = len(methods)
        heatmap_data = ones((num_methods, num_methods))

        # I know this is inefficient, but it really doesn't matter for what
        # we're doing here.
        for method1_idx, method1 in enumerate(methods):
            for method2_idx, method2 in enumerate(methods):
                corr_coeff = correlation_fn(method_data[method1],
                                            method_data[method2])
                heatmap_data[method1_idx][method2_idx] = corr_coeff

        heatmaps[correlation_name] = heatmap_data

    return heatmaps
