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
