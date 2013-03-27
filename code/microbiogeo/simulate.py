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
from os.path import basename, join
from qiime.filter import (filter_mapping_file_from_mapping_f,
                          filter_samples_from_otu_table)
from qiime.parse import parse_mapping_file_to_dict
from qiime.util import create_dir, MetadataMap
from random import randint, sample

from microbiogeo.parse import parse_mantel_results, parse_morans_i_results
from microbiogeo.util import run_command, run_parallel_jobs

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
    create_dir(out_dir)
    otu_table_fp = join(in_dir, tests['study'], 'otu_table.biom')
    map_fp = join(in_dir, tests['study'], 'map.txt')
    map_f = open(map_fp, 'U')
    category = tests['category'][0]

    otu_table_f = open(otu_table_fp, 'U')
    otu_table = parse_biom_table(otu_table_f)
    otu_table_f.seek(0)
    num_samps = len(otu_table.SampleIds)

    cmds = []
    for samp_size in tests['sample_sizes']:
        samp_size_dir = join(out_dir, '%d' % samp_size)
        create_dir(samp_size_dir)

        if samp_size <= num_samps:
            simsam_rep_num = 1

            subset_otu_table, subset_map_str = choose_gradient_subsets(
                    otu_table_f, map_f, category, samp_size)
            otu_table_f.seek(0)
            map_f.seek(0)

            subset_otu_table_fp = join(samp_size_dir, basename(otu_table_fp))
            subset_otu_table_f = open(subset_otu_table_fp, 'w')
            subset_otu_table.getBiomFormatJsonString('microbiogeo workflow',
                                                     subset_otu_table_f)
            subset_otu_table_f.close()

            subset_map_fp = join(samp_size_dir, basename(map_fp))
            subset_map_f = open(subset_map_fp, 'w')
            subset_map_f.write(subset_map_str)
            subset_map_f.close()

            for d in tests['dissim']:
                dissim_dir = join(samp_size_dir, '%r' % d)
                cmd = 'simsam.py -i %s -t %s -o %r -d %r -n %d -m %s;' % (
                        subset_otu_table_fp, tree_fp, dissim_dir, d, 1,
                        subset_map_fp)
                cmd += 'distance_matrix_from_mapping.py -i %s -c %s -o %s' % (
                        join(dissim_dir, 'map_n%d_d%r.txt' %
                        (simsam_rep_num, d)), category,
                        join(dissim_dir, '%s_dm.txt' % category))
                cmds.append(cmd)
        else:
            #simsam_rep_num = int(ceil(samp_size / num_samps))
            #otu_table = choose_gradient_subsets(otu_table_f, map_f, category,
            #                                    samp_size)
            # Make sure to choose subset afterwards.
            # TODO finish
            pass

    run_parallel_jobs(cmds, run_command)

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
        'category': ['Gradient', 'b', 'Gradient'],
        'methods': {
            'mantel': parse_mantel_results,
            'morans_i': parse_morans_i_results
        }
    }

    generate_gradient_simulated_data(in_dir, out_dir, gradient_tests, tree_fp)


if __name__ == "__main__":
    main()
