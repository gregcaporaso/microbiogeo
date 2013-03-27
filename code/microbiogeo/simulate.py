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

def choose_gradient_subsets(otu_table_f, map_f, category, num_total_samples):
    otu_table = parse_biom_table(otu_table_f)
    mdm, _ = parse_mapping_file_to_dict(map_f)

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

def process_gradient_simulated_data(in_dir, out_dir, methods, category,
                                    metric, num_perms):
    dm_fps = sorted(glob(join(in_dir, '%s_dm_*.txt' % metric)))
    map_fps = sorted(glob(join(in_dir, 'map_*.txt')))
    grad_dm_fps = sorted(glob(join(in_dir, '%s_dm_*.txt' % category)))
    assert len(dm_fps) == len(map_fps) and len(map_fps) == len(grad_dm_fps)

    cmds = []
    for method in methods:
        method_dir = join(out_dir, method)
        create_dir(method_dir)

        for dm_fp, map_fp, grad_dm_fp in zip(dm_fps, map_fps, grad_dm_fps):
            n, d = dm_fp.split('_dm_', 2)[1].split('.txt', 2)[0].split('_')
            n = int(n.split('n', 2)[1])
            d = float(d.split('d', 2)[1])
            results_dir = join(method_dir, 'n%d_d%f' % (n, d))

            if not has_results(results_dir):
                if method == 'mantel' or method == 'mantel_corr':
                    in_dm_fps = ','.join((dm_fp, grad_dm_fp))

                    cmds.append('compare_distance_matrices.py --method %s '
                                '-n %d -i %s -o %s' % (method, num_perms,
                                                       in_dm_fps, results_dir))
                elif method == 'morans_i':
                    cmds.append('compare_categories.py --method %s -i %s '
                                '-m %s -c %s -o %s' % (method, dm_fp, map_fp,
                                                       category, results_dir))

    run_parallel_jobs(cmds, run_command)

def create_sample_size_plots(in_dir, methods, category, metric, num_perms):
    for method, parse_fn in methods.items():
        method_dir = join(in_dir, method)

        results = defaultdict(dict)
        for res_dir in sorted(listdir(method_dir)):
            n, d = res_dir.split('_', 2)
            n = int(n.split('n', 2)[1])
            d = float(n.split('d', 2)[1])

def generate_gradient_simulated_data(in_dir, out_dir, tests, tree_fp):
    pass

def main():
    in_dir = 'test_datasets'
    out_dir = 'test_simulated_output'
    tree_fp = join('test_datasets', 'overview', 'rep_set.tre')
    gradient_tests = {
        'study': 'overview',
        'depth': 146,
        'metric': 'unweighted_unifrac',
        'num_perms': 999,
        'dissim': [0.001, 0.01, 0.1],
        'sample_sizes': [2, 3, 4, 5],
        'pos_control': ['Gradient', 'b', 'Gradient (positive control)'],
        'neg_control': ['DOB', 'r', 'Date of birth (negative control)'],
        'methods': {
            'mantel': parse_mantel_results,
            'morans_i': parse_morans_i_results
        }
    }

    generate_gradient_simulated_data(in_dir, out_dir, gradient_tests, tree_fp)


if __name__ == "__main__":
    main()
