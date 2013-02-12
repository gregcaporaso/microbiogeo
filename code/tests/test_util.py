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

from cogent.util.unit_test import TestCase, main
from qiime.parse import parse_distmat

from microbiogeo.util import (ExternalCommandFailedError, run_command,
                              shuffle_dm, subset_dm, subset_groups)

class UtilTests(TestCase):
    """Tests for the util.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.dm_f1 = dm_str1.split('\n')
        self.map_f1 = map_str1.split('\n')

    def test_run_command(self):
        """Test running an invalid command."""
        self.assertRaises(ExternalCommandFailedError, run_command,
                          'foobarbazbazbarfoo')

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


dm_str1 = """\tS1\tS2\tS3
S1\t0\t0.5\t0.7
S2\t0.5\t0\t0.1
S3\t0.7\t0.1\t0"""

map_str1 = """#SampleID\tBarcodeSequence\tCategory
S1\tAGCACGAGCCTA\tCat1
S2\tAGCACGAGCCTG\tCat2
S3\tAGCACGAGCCTC\tCat1
S4\tAGCACGAGCCTT\tCat1"""


if __name__ == "__main__":
    main()
