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

from microbiogeo.util import shuffle_dm, subsample_dm, subset_dm

class UtilTests(TestCase):
    """Tests for the util.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.dm_f1 = dm_str1.split('\n')

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


dm_str1 = """\tS1\tS2
S1\t0\t0.5
S2\t0.5\t0"""


if __name__ == "__main__":
    main()
