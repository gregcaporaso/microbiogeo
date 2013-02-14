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

from collections import defaultdict
from os import listdir
from os.path import exists
from random import sample, shuffle

from qiime.filter import filter_samples_from_distance_matrix
from qiime.format import format_distance_matrix
from qiime.parse import parse_distmat
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

def has_results(results_dir):
    """Returns True if results_dir exists and is not empty, False otherwise."""
    return exists(results_dir) and len(listdir(results_dir)) > 0

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

def is_empty(category_results):
    return (len(category_results) == 0) or \
           category_results['full'].isEmpty() or \
           category_results['shuffled'].isEmpty() or \
           [e for e in category_results['subsampled'] if e.isEmpty()]


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
