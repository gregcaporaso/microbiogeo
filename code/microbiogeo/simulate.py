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
from qiime.filter import (filter_mapping_file_from_mapping_f,
                          filter_samples_from_otu_table)
from qiime.parse import parse_mapping_file_to_dict
from qiime.util import MetadataMap
from random import randint, sample

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

def choose_gradient_subsets(otu_table_f, map_f, gradient, num_samples):
    otu_table = parse_biom_table(otu_table_f)
    mdm, _ = parse_mapping_file_to_dict(map_f)

    # Only keep the sample IDs that are in both the mapping file and OTU table.
    # Sort the samples according to the gradient.
    samp_ids = [(samp_id, float(metadata[gradient]))
                for samp_id, metadata in mdm.items()
                if samp_id in otu_table.SampleIds]
    samp_ids.sort(key=lambda samp_id: samp_id[1])

    samp_ids_to_keep = [samp_id[0] for samp_id in
                        _choose_evenly_spaced_items(samp_ids, num_samples)]

    assert len(samp_ids_to_keep) == num_samples, \
           "%d != %d" % (len(samp_ids_to_keep), num_samples)

    return (filter_samples_from_otu_table(otu_table, samp_ids_to_keep, 0, inf),
            filter_mapping_file_from_mapping_f(map_f, samp_ids_to_keep))

def _choose_evenly_spaced_items(sequence, num_items):
    # Adapted from http://stackoverflow.com/a/9873935
    return [sequence[int(ceil(i * float(len(sequence)) / num_items))]
            for i in range(num_items)]
