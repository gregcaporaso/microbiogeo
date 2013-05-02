#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Test suite for the util.py module."""

from os import chdir, getcwd
from os.path import exists, join
from shutil import rmtree
from tempfile import mkdtemp

from cogent.util.misc import remove_files
from cogent.util.unit_test import TestCase, main
from qiime.parse import parse_distmat
from qiime.util import get_qiime_temp_dir

from microbiogeo.util import (choose_gradient_subsets,
                              ExternalCommandFailedError, get_color_pool,
                              get_simsam_rep_num, has_results, is_empty,
                              run_command, run_parallel_jobs, shuffle_dm,
                              StatsResults, subset_dm, subset_groups)

class UtilTests(TestCase):
    """Tests for the util.py module functions."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.dm_f1 = dm_str1.split('\n')
        self.map_f1 = map_str1.split('\n')

        self.dm_f2 = dm_str2.split('\n')
        self.map_f2 = map_str2.split('\n')

        empty = StatsResults()
        nonempty = StatsResults()
        nonempty.addResult(0.1, 0.001)

        self.cat_res1 = {
            'full': empty,
            'shuffled': nonempty,
            'subsampled': [nonempty, nonempty]
        }

        self.cat_res2 = {
            'full': nonempty,
            'shuffled': empty,
            'subsampled': [nonempty, nonempty]
        }

        self.cat_res3 = {
            'full': nonempty,
            'shuffled': nonempty,
            'subsampled': [nonempty, empty]
        }

        self.cat_res4 = {
            'full': nonempty,
            'shuffled': nonempty,
            'subsampled': [nonempty, nonempty]
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

        # Set up temporary directories to use with tests.
        self.input_dir = mkdtemp(dir=self.tmp_dir,
                                 prefix='%s_input_dir_' % self.prefix)
        self.dirs_to_remove.append(self.input_dir)

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

    def test_run_command(self):
        """Test running an invalid command."""
        self.assertRaises(ExternalCommandFailedError, run_command,
                          'foobarbazbazbarfoo')

    def test_run_parallel_jobs(self):
        """Test running jobs in parallel."""
        # Doesn't error out if no jobs are submitted, which can happend during
        # a rerun of the workflow.
        self.assertTrue(run_parallel_jobs([], int) is None)

    def test_has_results(self):
        """Test checking a directory for results."""
        # Dir that doesn't exist.
        obs = has_results('/foobarbazbazbarfoo1234567890')
        self.assertFalse(obs)

        # Dir that exists but is empty.
        obs = has_results(self.input_dir)
        self.assertFalse(obs)

        # Dir that exists, with no required files, but is empty.
        obs = has_results(self.input_dir, required_files=[])
        self.assertFalse(obs)

        # Dir that exists and is not empty.
        tmp_fp = join(self.input_dir, 'foo.txt')
        tmp_f = open(tmp_fp, 'w')
        tmp_f.write('foo')
        tmp_f.close()
        self.files_to_remove.append(tmp_fp)

        obs = has_results(self.input_dir)
        self.assertTrue(obs)

        # Dir that exists and is not empty, with no required files.
        obs = has_results(self.input_dir, required_files=[])
        self.assertTrue(obs)

        # Dir that exists, is not empty, and has the required file.
        obs = has_results(self.input_dir, required_files=['foo.txt'])
        self.assertTrue(obs)

        # Dir that exists and is not empty, but doesn't have required files.
        obs = has_results(self.input_dir,
                          required_files=['foo.txt', 'bar.txt', 'baz.txt'])
        self.assertFalse(obs)

    def test_get_simsam_rep_num(self):
        """Test getting number of necessary simsam reps."""
        obs = get_simsam_rep_num(42, 13)
        self.assertEqual(obs, 4)

        obs = get_simsam_rep_num(50, 10)
        self.assertEqual(obs, 5)

        self.assertRaises(ValueError, get_simsam_rep_num, 42, 42)

    def test_shuffle_dm(self):
        """Test shuffling labels of distance matrix."""
        exp_labels, exp_dm = parse_distmat(self.dm_f1)

        order_changed = False
        for i in range(20):
            obs_labels, obs_dm = parse_distmat(
                    shuffle_dm(self.dm_f1).split('\n'))
            self.assertFloatEqual(obs_dm, exp_dm)

            try:
                self.assertIsPermutation(obs_labels, exp_labels)
            except AssertionError:
                pass
            else:
                order_changed = True

        self.assertTrue(order_changed)

    def test_subset_dm(self):
        """Test picking a subset of a distance matrix."""
        # Don't actually subset.
        exp = parse_distmat(self.dm_f1)
        obs = parse_distmat(subset_dm(self.dm_f1, 3).split('\n'))
        self.assertFloatEqual(obs, exp)

        obs_labels, obs_dm = parse_distmat(
                subset_dm(self.dm_f1, 1).split('\n'))
        self.assertEqual(len(obs_labels), 1)
        self.assertTrue(obs_labels[0] in exp[0])

        obs_labels, obs_dm = parse_distmat(
                subset_dm(self.dm_f1, 2).split('\n'))
        self.assertEqual(len(obs_labels), 2)
        self.assertTrue(obs_labels[0] in exp[0])
        self.assertTrue(obs_labels[1] in exp[0])

        self.assertRaises(ValueError, subset_dm, self.dm_f1, 4)

    def test_subset_groups(self):
        """Test picking subsets of sample groups in distance matrix."""
        # Don't filter anything out.
        exp = parse_distmat(self.dm_f1)
        obs = parse_distmat(subset_groups(
                self.dm_f1, self.map_f1, 'Category', 2).split('\n'))
        self.assertFloatEqual(obs, exp)

        obs = parse_distmat(subset_groups(
                self.dm_f1, self.map_f1, 'Category', 3).split('\n'))
        self.assertFloatEqual(obs, exp)

        # Pick groups of size 1.
        obs_labels, obs_dm = parse_distmat(subset_groups(
                self.dm_f1, self.map_f1, 'Category', 1).split('\n'))
        self.assertTrue('S2' in obs_labels)

        # XOR: either S1 or S3 should be in obs_labels, but not both.
        self.assertTrue(('S1' in obs_labels) != ('S3' in obs_labels))

    def test_choose_gradient_subsets(self):
        """Test picking subsets of gradients from a distance matrix."""
        # TODO test with size 3 (throws error).
        obs = choose_gradient_subsets(self.dm_f2, self.map_f2, 'Gradient',
                                      [2, 1, 5], 2)
        self.assertEqual(len(obs), 6)

        self.assertEqual(len(obs[0]), 2)
        self.assertEqual(len(obs[1]), 2)
        self.assertEqual(len(obs[2]), 1)
        self.assertEqual(len(obs[3]), 1)
        self.assertEqual(len(obs[4]), 5)
        self.assertEqual(len(obs[5]), 5)

        self.assertEqual(len(obs[0]), len(set(obs[0])))
        self.assertEqual(len(obs[1]), len(set(obs[1])))
        self.assertEqual(len(obs[2]), len(set(obs[2])))
        self.assertEqual(len(obs[3]), len(set(obs[3])))
        self.assertEqual(len(obs[4]), len(set(obs[4])))
        self.assertEqual(len(obs[5]), len(set(obs[5])))

        self.assertTrue(obs[0][0] in ['S5', 'S2'])
        self.assertTrue(obs[0][1] in ['S3', 'S1', 'S4'])

        self.assertTrue(obs[1][0] in ['S5', 'S2'])
        self.assertTrue(obs[1][1] in ['S3', 'S1', 'S4'])

        exp = ['S5', 'S2', 'S3', 'S1']
        self.assertTrue(obs[2][0] in exp)
        self.assertTrue(obs[3][0] in exp)

        exp = ['S5', 'S2', 'S3', 'S1', 'S4']
        self.assertEqual(obs[4], exp)
        self.assertEqual(obs[5], exp)

    def test_is_empty(self):
        """Test checking if category results are empty or not."""
        self.assertTrue(is_empty(self.cat_res1))
        self.assertTrue(is_empty(self.cat_res2))
        self.assertTrue(is_empty(self.cat_res3))
        self.assertTrue(is_empty({}))
        self.assertFalse(is_empty(self.cat_res4))

    def test_get_color_pool(self):
        """Test grabbing list of good colors to use."""
        obs = get_color_pool()
        self.assertEqual(len(obs), 27)

        obs2 = get_color_pool()
        self.assertFloatEqual(obs, obs2)
        self.assertFalse(obs is obs2)


class StatsResultsTests(TestCase):
    """Tests for the util.StatsResults class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.sr1 = StatsResults()

    def test_addResult(self):
        """Adding effect size and p-value works correctly on valid input."""
        self.sr1.addResult(0.5, 0.01)
        self.sr1.addResult(0.5, 0.001)
        self.assertFloatEqual(self.sr1.effect_size, 0.5)
        self.assertFloatEqual(self.sr1.p_values, [0.01, 0.001])

    def test_addResult_invalid_input(self):
        """Adding invalid input raises error."""
        # Effect sizes don't match.
        self.sr1.addResult(0.5, 0.01)
        self.assertRaises(ValueError, self.sr1.addResult, 0.6, 0.001)
        self.assertFloatEqual(self.sr1.effect_size, 0.5)
        self.assertFloatEqual(self.sr1.p_values, [0.01])

        # Invalid p-value range.
        self.sr1 = StatsResults()
        self.assertRaises(ValueError, self.sr1.addResult, 0.5, 1.1)
        self.assertTrue(self.sr1.effect_size is None)
        self.assertEqual(self.sr1.p_values, [])

        self.sr1.addResult(0.5, 0.01)
        self.sr1.addResult(0.5, 0.02)
        self.assertRaises(ValueError, self.sr1.addResult, 0.5, 1.1)
        self.assertRaises(ValueError, self.sr1.addResult, 0.5, -0.2)
        self.assertFloatEqual(self.sr1.effect_size, 0.5)
        self.assertFloatEqual(self.sr1.p_values, [0.01, 0.02])

    def test_isEmpty(self):
        """Test checking if results are empty or not."""
        self.assertTrue(self.sr1.isEmpty())

        self.sr1.addResult(0.5, 0.01)
        self.assertFalse(self.sr1.isEmpty())

    def test_str(self):
        """Test __str__ method."""
        # Empty results.
        obs = str(self.sr1)
        self.assertEqual(obs, 'Empty results')

        # Populated results.
        self.sr1.addResult(0.5, 0.01)
        self.sr1.addResult(0.5, 0.05)
        obs = str(self.sr1)
        self.assertEqual(obs, '0.50; ***, **')

    def test_check_p_value(self):
        """Raises error on invalid p-value."""
        self.sr1._check_p_value(0.0)
        self.sr1._check_p_value(0.5)
        self.sr1._check_p_value(1.0)

        self.assertRaises(ValueError, self.sr1._check_p_value, 1.5)
        self.assertRaises(ValueError, self.sr1._check_p_value, -1.5)

    def test_format_p_value_as_asterisk(self):
        """Test formatting a p-value to indicate statistical significance."""
        obs = self.sr1._format_p_value_as_asterisk(1.0)
        self.assertEqual(obs, 'x')

        obs = self.sr1._format_p_value_as_asterisk(0.09)
        self.assertEqual(obs, '*')

        obs = self.sr1._format_p_value_as_asterisk(0.045)
        self.assertEqual(obs, '**')

        obs = self.sr1._format_p_value_as_asterisk(0.01)
        self.assertEqual(obs, '***')

        obs = self.sr1._format_p_value_as_asterisk(0.0005)
        self.assertEqual(obs, '****')

    def test_format_p_value_as_asterisk_invalid_input(self):
        """Test supplying an invalid p-value results in error being thrown."""
        self.assertRaises(TypeError, self.sr1._format_p_value_as_asterisk, 1)
        self.assertRaises(TypeError, self.sr1._format_p_value_as_asterisk,
                          "0.05")
        self.assertRaises(TypeError, self.sr1._format_p_value_as_asterisk,
                          [0.05])

        self.assertRaises(ValueError, self.sr1._format_p_value_as_asterisk,
                          1.1)
        self.assertRaises(ValueError, self.sr1._format_p_value_as_asterisk,
                          -0.042)


dm_str1 = """\tS1\tS2\tS3
S1\t0\t0.5\t0.7
S2\t0.5\t0\t0.1
S3\t0.7\t0.1\t0"""

map_str1 = """#SampleID\tBarcodeSequence\tCategory\tGradient
S1\tAGCACGAGCCTA\tCat1\t4
S2\tAGCACGAGCCTG\tCat2\t2
S3\tAGCACGAGCCTC\tCat1\t3
S4\tAGCACGAGCCTT\tCat1\t2"""

dm_str2 = """\tS1\tS2\tS3\tS4\tS5
S1\t0\t0.5\t0.7\t0.2\t0.2
S2\t0.5\t0\t0.1\t0.4\t0.75
S3\t0.7\t0.1\t0\t0.1\t0.23
S4\t0.2\t0.4\t0.1\t0\t0.02
S5\t0.2\t0.75\t0.23\t0.02\t0"""

map_str2 = """#SampleID\tBarcodeSequence\tGradient
S1\tAGCACGAGCCTA\t4
S2\tAGCACGAGCCTG\t2
S3\tAGCACGAGCCTC\t3
S4\tAGCACGAGCCTT\t5
S5\tAGCACGAGCCTT\t1
S6\tAGCACGAGCCTT\t2"""


if __name__ == "__main__":
    main()
