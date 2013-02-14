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
from shutil import move

from qiime.util import create_dir

from microbiogeo.util import run_command, shuffle_dm, subset_dm, subset_groups

# The parameters need to be wrapped in parens in order to work with map. We put
# these functions in their own module to avoid import issues when using
# IPython.parallel. See http://stackoverflow.com/a/12307741 for more details.
def generate_per_study_depth_dms((in_dir, out_dir, study, depth, metrics,
                                  grouping_categories, gradient_categories,
                                  group_sizes, subset_sizes, num_shuffled,
                                  num_subsets, tree_fp)):
    in_study_dir = join(in_dir, study)
    out_study_dir = join(out_dir, study)
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
            shuffled_dm_f.write(shuffle_dm(renamed_dm_f))
            shuffled_dm_f.close()
            renamed_dm_f.close()

        # Create distance matrices with subsets of groups.
        for category in grouping_categories:
            for group_size in group_sizes:
                for i in range(1, num_subsets + 1):
                    subset_dm_fp = join(bdiv_out_dir,
                            '%s_dm_%s_gs%d_%d.txt' % (metric, category,
                                                      group_size, i))
                    subset_dm_f = open(subset_dm_fp, 'w')
                    subset_dm_f.write(subset_groups(
                            open(renamed_dm_fp, 'U'), open(map_fp, 'U'),
                            category, group_size))
                    subset_dm_f.close()

        # Filter the keyboard study distance matrix to include only samples
        # taken from keys of subjects M2, M3, and M9 before creating subsets.
        # We only want to include these samples (not human subject fingertips)
        # because we want to see if keys that are closer to each other are
        # correlated with community similarity.
        if study == 'keyboard':
            filtered_dm_fp = join(bdiv_out_dir,
                    add_filename_suffix(renamed_dm_fp, '_keys_only'))
            cmd = ('filter_distance_matrix.py -i %s, -o %s -m %s -s '
                   '\'COMMON_NAME:keyboard;HOST_SUBJECT_ID:M2,M3,M9\'' % (
                       renamed_dm_fp, filtered_dm_fp, map_fp))
            run_command(cmd)
            renamed_dm_fp = filtered_dm_fp

        # Create subsets of each non-shuffled distance matrix.
        for subset_size in subset_sizes:
            for i in range(1, num_subsets + 1):
                subset_dm_fp = join(bdiv_out_dir,
                        '%s_dm_n%d_%d.txt' % (metric, subset_size, i))
                subset_dm_f = open(subset_dm_fp, 'w')
                subset_dm_f.write(
                        subset_dm(open(renamed_dm_fp, 'U'), subset_size))
                subset_dm_f.close()
