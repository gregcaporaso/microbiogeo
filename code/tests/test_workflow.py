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

from microbiogeo.method import Adonis, Anosim, Mantel, MantelCorrelogram, Best
from microbiogeo.util import StatsResults
from microbiogeo.workflow import (_build_per_metric_real_data_commands,
                                  _collate_real_data_results,
                                  _parse_original_results_file,
                                  _parse_shuffled_results_files)

class WorkflowTests(TestCase):
    """Tests for the workflow.py module functions."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.cluster_workflow = {
            'overview': {
                'categories': [('Treatment', 'Treatment Category',
                                {'Control': 'Control', 'Fast': 'Fast'})],
                'depths': [(146, '25_percent')],
                'metrics': [('unweighted_unifrac', 'Unweighted UniFrac')],
                'num_real_data_perms': [99, 999],
                'num_sim_data_perms': 999,
                'dissim': [0.0, 0.001, 0.01, 0.1, 1.0, 10.0],
                'pcoa_dissim': [0.0, 0.001, 1.0, 10.0],
                'sample_sizes': [3, 5, 13],
                'pcoa_sample_size': 13,
                'num_sim_data_trials': 3,
                'num_shuffled_trials': 2,
                'methods': [Adonis(), Anosim()]
            }
        }

        self.gradient_workflow = {
            'overview': {
                'categories': [('DOB', 'Date of birth')],
                'depths': [(50, '2_percent'), (100, '5_percent'),
                           (146, '25_percent')],
                'metrics': [('unweighted_unifrac', 'Unweighted UniFrac'),
                            ('weighted_unifrac', 'Weighted UniFrac')],
                'num_real_data_perms': [99, 999],
                'num_sim_data_perms': 999,
                'dissim': [0.0, 0.001, 0.01, 0.1, 1.0, 10.0],
                'pcoa_dissim': [0.0, 0.001, 1.0, 10.0],
                'sample_sizes': [3, 5, 13],
                'pcoa_sample_size': 13,
                'num_sim_data_trials': 3,
                'num_shuffled_trials': 2,
                # Should both get skipped during collation.
                'methods': [MantelCorrelogram(), Best()],
                'best_method_env_vars': ['DOB']
            }
        }

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

    def test_collate_real_data_results(self):
        """Test collating method results."""
        # These methods should be skipped.
        obs = _collate_real_data_results('/foobarbaz123',
                                         self.gradient_workflow)
        self.assertEqual(obs, exp_collate_real_data_results1)

        # Bogus paths/input, so should get empty results back.
        obs = _collate_real_data_results('/foobarbaz123',
                                         self.cluster_workflow)

        inner_obs = obs['25_percent']['unweighted_unifrac']['anosim']['overview']['Treatment']
        self.assertTrue(inner_obs['original'].isEmpty())
        self.assertTrue(inner_obs['shuffled'].isEmpty())

        inner_obs = obs['25_percent']['unweighted_unifrac']['adonis']['overview']['Treatment']
        self.assertTrue(inner_obs['original'].isEmpty())
        self.assertTrue(inner_obs['shuffled'].isEmpty())

    def test_parse_original_results_file(self):
        res = StatsResults()
        _parse_original_results_file('/foobarbaz123', Anosim(), 'Treatment',
                                     res)
        self.assertTrue(res.isEmpty())

        _parse_original_results_file('/foobarbaz123', Anosim(), 'Treatment',
                                     res, 42)
        self.assertTrue(res.isEmpty())

    def test_parse_shuffled_results_files(self):
        res = StatsResults()
        _parse_shuffled_results_files('/foobarbaz123', Adonis(), 'Treatment',
                                     res, 10)
        self.assertTrue(res.isEmpty())

        _parse_shuffled_results_files('/foobarbaz123', Anosim(), 'Treatment',
                                     res, 20, 88)
        self.assertTrue(res.isEmpty())


exp_collate_real_data_results1 = {'5_percent': {'weighted_unifrac': {}, 'unweighted_unifrac': {}}, '25_percent': {'weighted_unifrac': {}, 'unweighted_unifrac': {}}, '2_percent': {'weighted_unifrac': {}, 'unweighted_unifrac': {}}}


if __name__ == "__main__":
    main()
