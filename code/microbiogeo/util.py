#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Module with various utility functions."""

from biom.parse import parse_biom_table
from collections import defaultdict
from os import listdir
from os.path import exists, join
from random import randint, sample, shuffle

from IPython.parallel import Client

from numpy import ceil

from qiime.colors import data_colors, data_color_order
from qiime.filter import filter_samples_from_distance_matrix
from qiime.format import format_distance_matrix
from qiime.make_distance_histograms import matplotlib_rgb_color
from qiime.parse import parse_distmat, parse_mapping_file_to_dict
from qiime.util import MetadataMap, qiime_system_call

class ExternalCommandFailedError(Exception):
    pass

def run_command(cmd):
    stdout, stderr, ret_val = qiime_system_call(cmd)

    if ret_val != 0:
        raise ExternalCommandFailedError("The command '%s' failed with exit "
                                         "status %d.\n\nStdout:\n\n%s\n\n"
                                         "Stderr:\n\n%s\n" % (cmd,
                                         ret_val, stdout, stderr))

def run_parallel_jobs(jobs, job_fn):
    # IPython will error out if jobs is empty.
    if jobs:
        c = Client()
        lview = c.load_balanced_view()
        lview.block = True
        lview.map(job_fn, jobs)

def has_results(results_dir, required_files=None):
    """Returns True if results_dir exists and is not empty, False otherwise.
    
    If required_files is provided, each filename in the list must be present in
    results_dir in order for this function to return True. Otherwise, a simple
    check that the directory exists and isn't empty is performed.
    """
    has_results = exists(results_dir) and len(listdir(results_dir)) > 0

    if has_results and required_files:
        for req_file in required_files:
            if not exists(join(results_dir, req_file)):
                has_results = False
                break

    return has_results

def get_num_samples(table_fp):
    """Returns the number of samples in the table."""
    with open(table_fp, 'U') as table_f:
        table = parse_biom_table(table_f)
        return len(table.SampleIds)

def shuffle_dm(dm_f):
    labels, dm_data = parse_distmat(dm_f)
    shuffle(labels)
    return format_distance_matrix(labels, dm_data)

def subset_dm(dm_f, num_samps):
    labels, dm_data = parse_distmat(dm_f)
    samp_ids_to_keep = sample(labels, num_samps)
    return filter_samples_from_distance_matrix((labels, dm_data),
                                               samp_ids_to_keep, negate=True)

def subset_groups(dm_f, map_f, category, max_group_size):
    dm_labels, dm_data = parse_distmat(dm_f)
    metadata_map = MetadataMap.parseMetadataMap(map_f)

    category_map = defaultdict(list)
    for samp_id in metadata_map.SampleIds:
        # Mapping files can have more samples than distance matrices, which can
        # happen in this case since we are dealing with rarefied OTU tables
        # (samples get dropped).
        if samp_id in dm_labels:
            category_val = metadata_map.getCategoryValue(samp_id, category)
            category_map[category_val].append(samp_id)

    samp_ids_to_keep = []
    for category_val, samp_ids in category_map.items():
        samp_ids_to_keep.extend(
                sample(samp_ids, min(max_group_size, len(samp_ids))))

    return filter_samples_from_distance_matrix((dm_labels, dm_data),
                                               samp_ids_to_keep, negate=True)

def choose_gradient_subsets(dm_f, map_f, gradient, subset_sizes, num_subsets):
    subsets = []

    mdm, _ = parse_mapping_file_to_dict(map_f)
    dm_labels, dm_data = parse_distmat(dm_f)

    # Only keep the sample IDs that are in both the mapping file and distance
    # matrix.
    samp_ids = [(samp_id, float(metadata[gradient]))
                for samp_id, metadata in mdm.items() if samp_id in dm_labels]
    samp_ids.sort(key=lambda samp_id: samp_id[1])

    for subset_size in subset_sizes:
        # Adapted from http://stackoverflow.com/a/9873935
        # We add 1 to the number of samples we want because we want subset_size
        # intervals to choose from.
        bin_idxs = [int(ceil(i * len(samp_ids) / (subset_size + 1)))
                    for i in range(subset_size + 1)]

        for subset_num in range(num_subsets):
            samp_ids_to_keep = []

            for i in range(len(bin_idxs) - 1):
                if i == len(bin_idxs) - 2:
                    # We're at the last bin, so choose from the entire bin
                    # range.
                    if bin_idxs[i + 1] < len(samp_ids):
                        end_idx = bin_idxs[i + 1]
                    else:
                        end_idx = bin_idxs[i + 1] - 1

                    samp_ids_to_keep.append(
                            samp_ids[randint(bin_idxs[i], end_idx)][0])
                else:
                    # We subtract 1 since randint is inclusive on both sides,
                    # and we don't want to choose the same sample ID multiple
                    # times from different bins.
                    samp_ids_to_keep.append(
                            samp_ids[randint(bin_idxs[i],
                                             bin_idxs[i + 1] - 1)][0])

            assert len(samp_ids_to_keep) == subset_size, \
                   "%d != %d" % (len(samp_ids_to_keep), subset_size)

            subsets.append(samp_ids_to_keep)

    return subsets

def is_empty(category_results):
    return (len(category_results) == 0) or \
           category_results['full'].isEmpty() or \
           category_results['shuffled'].isEmpty() or \
           [e for e in category_results['subsampled'] if e.isEmpty()]

def get_color_pool():
    # We don't like yellow...
    color_order = data_color_order[:]
    color_order.remove('yellow1')
    color_order.remove('yellow2')
    return [matplotlib_rgb_color(data_colors[color].toRGB())
            for color in color_order]

class StatsResults(object):

    def __init__(self):
        self.effect_size = None
        self.p_values = []

    def addResult(self, effect_size, p_value):
        self._check_p_value(p_value)

        if self.isEmpty():
            self.effect_size = effect_size
        else:
            if effect_size != self.effect_size:
                raise ValueError("The effect size %.4f is not the same as the "
                                 "previously supplied effect size %.4f. These "
                                 "must be the same for different numbers of "
                                 "permutations." % (effect_size,
                                                    self.effect_size))

        self.p_values.append(p_value)

    def isEmpty(self):
        return self.effect_size is None

    def __str__(self):
        if self.isEmpty():
            result = 'Empty results'
        else:
            result = '%.2f; %s' % (self.effect_size, ', '.join(
                    map(self._format_p_value_as_asterisk, self.p_values)))
        return result

    def _check_p_value(self, p_value):
        if p_value < 0 or p_value > 1:
            raise ValueError("Invalid p-value: %.4f" % p_value)

    def _format_p_value_as_asterisk(self, p_value):
        if not isinstance(p_value, float):
            raise TypeError("p-value must be a float.")
        self._check_p_value(p_value)

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
