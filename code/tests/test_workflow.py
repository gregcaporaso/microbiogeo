#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Test suite for the workflow.py module."""

from os import chdir, getcwd
from os.path import exists, join
from shutil import rmtree
from tempfile import mkdtemp

from cogent.util.misc import remove_files
from cogent.util.unit_test import TestCase, main
from qiime.util import create_dir, get_qiime_temp_dir

from microbiogeo.method import Adonis, Mantel, MantelCorrelogram, Best
from microbiogeo.util import StatsResults
from microbiogeo.workflow import (_build_per_metric_real_data_commands,
                                  _collate_category_results, _collate_results)

class WorkflowTests(TestCase):
    """Tests for the workflow.py module functions."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.studies1 = {
            'overview': {
                'depths': [50, 100, 146],
                'grouping_categories': ['Treatment'],
                'gradient_categories': ['DOB'],
                'group_sizes': [3, 4],
                'subset_sizes': [3, 4],
                'best_method_env_vars': ['DOB']
            }
        }

        # Invalid method type.
        self.methods1 = {
            'foo': [Adonis()]
        }

        # Methods that should be skipped during collation.
        self.methods2 = {
            'gradient': [MantelCorrelogram(), Best()]
        }

        # The prefix to use for temporary files/dirs. This prefix may be added
        # to, but all temp dirs and files created by the tests will have this
        # prefix at a minimum.
        self.prefix = 'microbiogeo_tests'

        self.start_dir = getcwd()
        self.dirs_to_remove = []
        self.files_to_remove = []

        self.tmp_dir = get_qiime_temp_dir()

        if not exists(self.tmp_dir):
            makedirs(self.tmp_dir)

            # If test creates the temp dir, also remove it.
            self.dirs_to_remove.append(self.tmp_dir)

        # Set up a temporary input directory structure:
        #   <tmp location>/overview/bdiv_even50/dm.txt
        self.input_dir = mkdtemp(dir=self.tmp_dir,
                                 prefix='%s_input_dir_' % self.prefix)
        self.dirs_to_remove.append(self.input_dir)

        study_dir = join(self.input_dir, 'overview')
        create_dir(study_dir)
        self.dirs_to_remove.append(study_dir)

        depth_dir = join(study_dir, 'bdiv_even50')
        create_dir(depth_dir)
        self.dirs_to_remove.append(depth_dir)

        tmp_fp = join(depth_dir, 'bray_curtis_dm.txt')
        tmp_f = open(tmp_fp, 'w')
        tmp_f.write('foo')
        tmp_f.close()
        self.files_to_remove.append(tmp_fp)

    def tearDown(self):
        """Remove temporary files/dirs created by tests."""
        # Change back to the start dir - some workflows change directory.
        chdir(self.start_dir)
        remove_files(self.files_to_remove)

        # Remove directories last, so we don't get errors trying to remove
        # files which may be in the directories.
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

    def test_build_per_metric_real_data_commands(self):
        exp = 'beta_diversity.py -i /bar/baz.biom -o /foo/original -m unweighted_unifrac -t /bar/tree.tre && mv /foo/original/unweighted_unifrac_baz.txt /foo/original/dm.txt && cp /map.txt /foo/original/map.txt && principal_coordinates.py -i /foo/original/dm.txt -o /foo/original/pc.txt && mkdir -p /foo/0 && shuffle_distance_matrix.py -i /foo/original/dm.txt -o /foo/0/dm.txt && cp /foo/original/map.txt /foo/0/map.txt && principal_coordinates.py -i /foo/0/dm.txt -o /foo/0/pc.txt && mkdir -p /foo/1 && shuffle_distance_matrix.py -i /foo/original/dm.txt -o /foo/1/dm.txt && cp /foo/original/map.txt /foo/1/map.txt && principal_coordinates.py -i /foo/1/dm.txt -o /foo/1/pc.txt'
        obs = _build_per_metric_real_data_commands('cluster', '/foo',
                '/bar/baz.biom', '/map.txt', '/bar/tree.tre',
                ('unweighted_unifrac', 'Unweighted UniFrac'), ['A', 'B'], 2)
        self.assertEqual(obs, exp)

        exp = 'beta_diversity.py -i /bar/baz.biom -o /foo/original -m unweighted_unifrac -t /bar/tree.tre && mv /foo/original/unweighted_unifrac_baz.txt /foo/original/dm.txt && cp /map.txt /foo/original/map.txt && principal_coordinates.py -i /foo/original/dm.txt -o /foo/original/pc.txt && distance_matrix_from_mapping.py -i /foo/original/map.txt -c A -o /foo/original/A_dm.txt && distance_matrix_from_mapping.py -i /foo/original/map.txt -c B -o /foo/original/B_dm.txt && mkdir -p /foo/0 && shuffle_distance_matrix.py -i /foo/original/dm.txt -o /foo/0/dm.txt && cp /foo/original/map.txt /foo/0/map.txt && principal_coordinates.py -i /foo/0/dm.txt -o /foo/0/pc.txt && distance_matrix_from_mapping.py -i /foo/0/map.txt -c A -o /foo/0/A_dm.txt && distance_matrix_from_mapping.py -i /foo/0/map.txt -c B -o /foo/0/B_dm.txt && mkdir -p /foo/1 && shuffle_distance_matrix.py -i /foo/original/dm.txt -o /foo/1/dm.txt && cp /foo/original/map.txt /foo/1/map.txt && principal_coordinates.py -i /foo/1/dm.txt -o /foo/1/pc.txt && distance_matrix_from_mapping.py -i /foo/1/map.txt -c A -o /foo/1/A_dm.txt && distance_matrix_from_mapping.py -i /foo/1/map.txt -c B -o /foo/1/B_dm.txt'
        obs = _build_per_metric_real_data_commands('gradient', '/foo',
                '/bar/baz.biom', '/map.txt', '/bar/tree.tre',
                ('unweighted_unifrac', 'Unweighted UniFrac'), ['A', 'B'], 2)
        self.assertEqual(obs, exp)

    def test_collate_results(self):
        """Test collating method results."""
        # These methods should be skipped.
        obs = _collate_results(self.input_dir, self.studies1, self.methods2,
                ['5_percent', '25_percent'], ['euclidean', 'bray_curtis'],
                [42, 43], 2, 3)
        self.assertEqual(obs, exp_collate_results1)

    def test_collate_results_invalid_input(self):
        """Test collating method results with invalid input raises an error."""
        # Invalid method type.
        self.assertRaises(ValueError, _collate_results, self.input_dir,
                self.studies1, self.methods1, ['5_percent', '25_percent'],
                ['euclidean', 'bray_curtis'], [42, 43], 2, 3)

    def test_collate_category_results(self):
        """Test collating category results."""
        # Bogus paths/input, so should get empty results back.
        full_results = StatsResults()
        shuffled_results = StatsResults()
        ss_results = [StatsResults()]

        _collate_category_results(full_results, shuffled_results, ss_results,
                '/foobarbaz123', 'foo', 42, 'euclidean', 'gradient', Mantel(),
                'DOB', [2], 3, 4, permutation=47)

        self.assertTrue(full_results.isEmpty())
        self.assertTrue(shuffled_results.isEmpty())
        self.assertEqual(len(ss_results), 1)
        self.assertTrue(ss_results[0].isEmpty())

    def test_collate_category_results_invalid_input(self):
        """Test collating category results w/ invalid input raises an error."""
        # Invalid method type.
        self.assertRaises(ValueError, _collate_category_results,
                StatsResults(), StatsResults(), [StatsResults()],
                '/foobarbaz123', 'foo', 42, 'euclidean', 'foo', Mantel(),
                'DOB', [2], 3, 4, permutation=47)


exp_collate_results1 = {'5_percent': {'euclidean': {'gradient': {}}, 'bray_curtis': {'gradient': {}}}, '25_percent': {'euclidean': {'gradient': {}}, 'bray_curtis': {'gradient': {}}}}


if __name__ == "__main__":
    main()
