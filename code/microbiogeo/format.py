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

def format_method_comparison_table(methods_results):
    constructed_header = None
    rows = []

    for method, method_res in sorted(methods_results.items()):
        header = ['Method']
        row = [method]

        for study, study_res in sorted(method_res.items()):
            for category, category_res in sorted(study_res.items()):
                header_column_title = '%s\r%s' % (study, category)
                header.append(header_column_title)
                header.append(header_column_title + ' (shuffled)')

                if len(category_res) == 0:
                    row.extend(['N/A'] * 2)
                else:
                    if category_res['original'].isEmpty():
                        row.append('N/A')
                    else:
                        row.append(str(category_res['original']))

                    if category_res['shuffled'].isEmpty():
                        row.append('N/A')
                    else:
                        row.append(str(category_res['shuffled']))

        if constructed_header is None:
            rows.append(header)
            constructed_header = header
        elif constructed_header != header:
            raise ValueError("The studies and/or categories did not match up "
                             "exactly between one or more of the methods.")

        rows.append(row)

    return rows

def format_method_comparison_heatmaps(real_data_results, sim_data_results,
                                      heatmap_methods):
    shared_studies = None
    shared_categories = {}

    for depth_desc, depth_res in real_data_results.items():
        for metric, metric_res in depth_res.items():
            for method, method_res in metric_res.items():
                matched_method = False
                for m in heatmap_methods:
                    if method == m.Name:
                        matched_method = True
                        break

                if matched_method:
                    studies = sorted(method_res.keys())

                    if shared_studies is None:
                        shared_studies = studies
                    elif studies != shared_studies:
                        raise ValueError("Not all methods to include in "
                                         "the heatmap have results for "
                                         "the same studies.")

                    for study, study_res in sorted(method_res.items()):
                        categories = [cat for cat, cat_res in \
                                      sorted(study_res.items()) if not
                                      is_empty(cat_res)]

                        if study not in shared_categories:
                            shared_categories[study] = set(categories)
                        else:
                            shared_categories[study] &= set(categories)

    # Gather real data effect sizes for each method (in the same order for each
    # method!).
    method_data = defaultdict(list)
    for depth_desc, depth_res in real_data_results.items():
        for metric, metric_res in depth_res.items():
            for method, method_res in metric_res.items():
                matched_method = False
                for m in heatmap_methods:
                    if method == m.Name:
                        matched_method = True
                        break

                if not matched_method:
                    continue

                for study, study_res in sorted(method_res.items()):
                    for category, category_res in sorted(study_res.items()):
                        if category in shared_categories[study]:
                            method_data[method].append(
                                    category_res['original'].effect_size)
                            method_data[method].append(
                                    category_res['shuffled'].effect_size)

    # Gather simulated data effect sizes.
    for method, method_res in sim_data_results.items():
        matched_method = False
        for m in heatmap_methods:
            if method == m.Name:
                matched_method = True
                break

        if not matched_method:
            continue

        for study, study_res in sorted(method_res.items()):
            for depth, depth_res in sorted(study_res.items()):
                for category, category_res in sorted(depth_res.items()):
                    for trial_num, trial_num_res in sorted(category_res.items()):
                        for samp_size, samp_size_res in sorted(trial_num_res.items()):
                            for dissim, dissim_res in sorted(samp_size_res.items()):
                                for metric, metric_res in sorted(dissim_res.items()):
                                    if metric_res.isEmpty():
                                        raise ValueError("Encountered empty "
                                                "simulated data results.")
                                    else:
                                        method_data[method].append(
                                                metric_res.effect_size)

    # Make sure our data looks sane. We should have the same number of
    # observations (i.e. effect sizes) for each method.
    data_length = None
    for method, data in method_data.items():
        if data_length is None:
            data_length = len(data)
        elif len(data) != data_length:
            raise ValueError("The number of observations (i.e. effect sizes) "
                             "is not the same between all methods, so we "
                             "can't compare them.")

    # Compute the correlation coefficient between each pair of methods and put
    # the output in an array. This array can then be used to generate a
    # text-based table or heatmap.
    results = {}
    for correlation_name, correlation_fn in \
            ('pearson', pearson), ('spearman', spearman):
        num_methods = len(heatmap_methods)
        heatmap_data = ones((num_methods, num_methods))

        # I know this is inefficient, but it really doesn't matter for what
        # we're doing here.
        for method1_idx, method1 in enumerate(heatmap_methods):
            for method2_idx, method2 in enumerate(heatmap_methods):
                corr_coeff = correlation_fn(method_data[method1.Name],
                                            method_data[method2.Name])
                heatmap_data[method1_idx][method2_idx] = corr_coeff

        # Mask out the upper triangle. Taken from
        # http://stackoverflow.com/a/2332520
        mask = invert(tri(heatmap_data.shape[0], k=0, dtype=bool))
        heatmap_data = numpy.ma.array(heatmap_data, mask=mask)

        results[correlation_name] = heatmap_data

    return results
