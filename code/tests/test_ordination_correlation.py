#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Test suite for the ordination_correlation.py module."""

from cogent.util.unit_test import TestCase, main

from microbiogeo.ordination_correlation import compute_ordination_correlation

class OrdinationCorrelationTests(TestCase):
    """Tests for the ordination_correlation.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.coords1 = coords1.split('\n')
        self.map1 = map1.split('\n')
        self.category1 = 'Gradient'

    def test_compute_ordination_correlation(self):
        """Functions correctly using valid input."""
        # Defaults.
        exp = (1.0, 0.0)
        obs = compute_ordination_correlation(self.map1, self.coords1,
                                             self.category1)
        self.assertFloatEqual(obs[:2], exp)
        self.assertIsProb(float(obs[2]))

        # Spearman with second axis.
        exp = (0.4, 0.6)
        obs = compute_ordination_correlation(self.map1, self.coords1,
                                             self.category1,
                                             axis=2,
                                             correlation_type='spearman')
        self.assertFloatEqual(obs[:2], exp)
        self.assertIsProb(float(obs[2]))

        # No permutations.
        exp = (1.0, 0.0, 'N/A')
        obs = compute_ordination_correlation(self.map1, self.coords1,
                                             self.category1,
                                             num_permutations=0)
        self.assertFloatEqual(obs, exp)

        # Few permutations.
        exp = (1.0, 0.0, 'Too few iters to compute p-value (num_iters=1)')
        obs = compute_ordination_correlation(self.map1, self.coords1,
                                             self.category1,
                                             num_permutations=1)
        self.assertEqual(obs, exp)

    def test_compute_ordination_correlation_invalid_input(self):
        """Correctly handles invalid input."""
        # category
        with self.assertRaises(ValueError):
            compute_ordination_correlation(self.map1, self.coords1, 'foo')

        # category state
        with self.assertRaises(ValueError):
            compute_ordination_correlation(self.map1, self.coords1,
                                           'Treatment')

        # axis
        with self.assertRaises(ValueError):
            compute_ordination_correlation(self.map1, self.coords1,
                                           self.category1, axis=-1)
        with self.assertRaises(ValueError):
            compute_ordination_correlation(self.map1, self.coords1,
                                           self.category1, axis=4)

        # correlation_type
        with self.assertRaises(ValueError):
            compute_ordination_correlation(self.map1, self.coords1,
                                           self.category1,
                                           correlation_type='foo')

        # num_permutations
        with self.assertRaises(ValueError):
            compute_ordination_correlation(self.map1, self.coords1,
                                           self.category1,
                                           num_permutations=-1)


coords1 = """pc vector number\t1\t2\t3
S3\t1.25\t0.53\t0.01
S2\t1.0\t0.03\t0.17
S1\t1.5\t0.08\t0.07
S4\t2.0\t0.12\t0.42


eigvals\t0.275632711048\t0.0786886203713\t0.0369532412336
% variation explained\t60.9625169162\t17.4037991799\t8.17305966788
"""

map1 = """#SampleID\tGradient\tTreatment
S1\t4.5\tA
S2\t4.0\tB
S3\t4.25\tA
S4\t5.0\tA
S5\t42.0\tB
"""


if __name__ == "__main__":
    main()
