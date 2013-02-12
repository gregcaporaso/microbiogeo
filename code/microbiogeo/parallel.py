#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Module for functionality meant to be run in parallel."""

from os.path import join
from shutil import copy, move

from qiime.util import create_dir

from microbiogeo.util import run_command

# The parameters need to be wrapped in parens in order to work with map.
def generate_per_study_depth_dms((study, depth, metrics, categories,
                                  group_sizes, num_shuffled, num_subsets,
                                  in_dir, out_dir, tree_fp, shuffle_dm_fn,
                                  subsample_dm_fn)):
    in_study_dir = join(in_dir, study)
    out_study_dir = join(out_dir, study)
    create_dir(out_study_dir)
    copy(join(in_study_dir, 'map.txt'), out_study_dir)
    map_fp = join(out_study_dir, 'map.txt')

    full_otu_fp = join(in_study_dir, 'otu_table.biom')
    even_otu_fp = join(out_study_dir, 'otu_table_even%d.biom' % depth)
    bdiv_out_dir = join(out_study_dir, 'bdiv_even%d' % depth)

    metrics_param = ','.join(metrics)

    cmd = 'single_rarefaction.py -i %s -o %s -d %d' % (full_otu_fp,
                                                       even_otu_fp,
                                                       depth)
    run_command(cmd)

    cmd = 'beta_diversity.py -i %s -o %s -m %s -t %s' % (even_otu_fp,
                                                         bdiv_out_dir,
                                                         metrics_param,
                                                         tree_fp)
    run_command(cmd)

    # Rename each file to match QIIME's standard naming conventions of
    # distance matrices. Generate shuffled versions of each distance matrix.
    for metric in metrics:
        dm_fp = join(bdiv_out_dir, '%s_otu_table_even%d.txt' % (metric, depth))
        renamed_dm_fp = join(bdiv_out_dir, '%s_dm.txt' % metric)
        move(dm_fp, renamed_dm_fp)

        for i in range(1, num_shuffled + 1):
            renamed_dm_f = open(renamed_dm_fp, 'U')
            shuffled_dm_fp = join(bdiv_out_dir,
                                  '%s_dm_shuffled%d.txt' % (metric, i))
            shuffled_dm_f = open(shuffled_dm_fp, 'w')
            shuffled_dm_f.write(shuffle_dm_fn(renamed_dm_f))
            shuffled_dm_f.close()
            renamed_dm_f.close()

        # Create subsampled distance matrices.
        for category in categories:
            for group_size in group_sizes:
                for i in range(1, num_subsets + 1):
                    subsampled_dm_fp = join(bdiv_out_dir,
                            '%s_dm_%s_ss%d_%d.txt' % (metric, category,
                                                      group_size, i))
                    subsampled_dm = open(subsampled_dm_fp, 'w')
                    subsampled_dm.write(subsample_dm_fn(
                            open(renamed_dm_fp, 'U'), open(map_fp, 'U'),
                            category, group_size))
                    subsampled_dm.close()
