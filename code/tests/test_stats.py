#!/usr/bin/env python
from __future__ import division

__author__ = "Michael Dwan"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Michael Dwan", "Logan Knecht", "Jai Ram Rideout",
               "Andrew Cochran"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Michael Dwan"
__email__ = "mdwan.tgen@gmail.com"
__status__ = "Development"

"""Test suite for classes, methods and functions of the stats module."""

from cogent.util.unit_test import TestCase, main
from numpy import array, matrix, roll, asarray
from numpy.random import permutation

from qiime.util import DistanceMatrix, MetadataMap
from python.qiime.stats import DistanceBasedRda

class TestHelper(TestCase):
    """Helper class that instantiates some commonly-used objects.

    This class should be subclassed by any test classes that want to use its
    members.
    """

    def setUp(self):
        """Define some useful test objects."""
        # The unweighted unifrac distance matrix from the overview tutorial.
        self.overview_dm_str = ["\tPC.354\tPC.355\tPC.356\tPC.481\tPC.593\
                                \tPC.607\tPC.634\tPC.635\tPC.636",
                                "PC.354\t0.0\t0.595483768391\t0.618074717633\
                                \t0.582763100909\t0.566949022108\
                                \t0.714717232268\t0.772001731764\
                                \t0.690237118413\t0.740681707488",
                                "PC.355\t0.595483768391\t0.0\t0.581427669668\
                                \t0.613726772383\t0.65945132763\
                                \t0.745176523638\t0.733836123821\
                                \t0.720305073505\t0.680785600439",
                                "PC.356\t0.618074717633\t0.581427669668\t0.0\
                                \t0.672149021573\t0.699416863323\
                                \t0.71405573754\t0.759178215168\
                                \t0.689701276341\t0.725100672826",
                                "PC.481\t0.582763100909\t0.613726772383\
                                \t0.672149021573\t0.0\t0.64756120797\
                                \t0.666018240373\t0.66532968784\
                                \t0.650464714994\t0.632524644216",
                                "PC.593\t0.566949022108\t0.65945132763\
                                \t0.699416863323\t0.64756120797\t0.0\
                                \t0.703720200713\t0.748240937349\
                                \t0.73416971958\t0.727154987937",
                                "PC.607\t0.714717232268\t0.745176523638\
                                \t0.71405573754\t0.666018240373\
                                \t0.703720200713\t0.0\t0.707316869557\
                                \t0.636288883818\t0.699880573956",
                                "PC.634\t0.772001731764\t0.733836123821\
                                \t0.759178215168\t0.66532968784\
                                \t0.748240937349\t0.707316869557\t0.0\
                                \t0.565875193399\t0.560605525642",
                                "PC.635\t0.690237118413\t0.720305073505\
                                \t0.689701276341\t0.650464714994\
                                \t0.73416971958\t0.636288883818\
                                \t0.565875193399\t0.0\t0.575788039321",
                                "PC.636\t0.740681707488\t0.680785600439\
                                \t0.725100672826\t0.632524644216\
                                \t0.727154987937\t0.699880573956\
                                \t0.560605525642\t0.575788039321\t0.0"]
        self.overview_dm = DistanceMatrix.parseDistanceMatrix(
            self.overview_dm_str)

        # The overview tutorial's metadata mapping file.
        self.overview_map_str = ["#SampleID\tBarcodeSequence\tTreatment\tDOB",
                                 "PC.354\tAGCACGAGCCTA\tControl\t20061218",
                                 "PC.355\tAACTCGTCGATG\tControl\t20061218",
                                 "PC.356\tACAGACCACTCA\tControl\t20061126",
                                 "PC.481\tACCAGCGACTAG\tControl\t20070314",
                                 "PC.593\tAGCAGCACTTGT\tControl\t20071210",
                                 "PC.607\tAACTGTGCGTAC\tFast\t20071112",
                                 "PC.634\tACAGAGTCGGCT\tFast\t20080116",
                                 "PC.635\tACCGCAGAGTCA\tFast\t20080116",
                                 "PC.636\tACGGTGAGTGTC\tFast\t20080116"]
        self.overview_map = MetadataMap.parseMetadataMap(self.overview_map_str)

        # A 1x1 dm.
        self.single_ele_dm = DistanceMatrix(array([[0]]), ['s1'], ['s1'])


class DistanceBasedRdaTests(TestHelper):
    """Tests for the DistanceBasedRda class."""

    def setUp(self):
        """Define some useful data to use in testing."""
        super(DistanceBasedRdaTests, self).setUp()
        self.dbrda = DistanceBasedRda(self.overview_dm,
                                      self.overview_map, "Treatment")

    def test_call(self):
        """Test running RDA over various inputs."""
        self.dbrda()

    def test_center_matrix(self):
        """Test the centering of matrices."""
        exp = matrix([[-0.5, -1.0, 2.5], [0.5, 1.0, -2.5]])
        obs = self.dbrda._center_matrix(matrix([[1, 2, 5], [2, 4, 0]]))
        self.assertFloatEqual(exp, obs)

    def test_create_factor(self):
        """Test creating factors from a group membership list."""
        cat_data = ["Fast", "Fast", "Control", "Fast", "Control"]
        exp = matrix([0, 0, 1, 0, 1]).T
        obs = self.dbrda._create_factor(cat_data)
        self.assertFloatEqual(exp, obs)

        num_data = [1, 2.0, 5, 0, -20, 99.99]
        exp = matrix(num_data).T
        obs = self.dbrda._create_factor(num_data)
        self.assertFloatEqual(exp, obs)

        small_data = ["foo"]
        exp = matrix([0]).T
        obs = self.dbrda._create_factor(small_data)
        self.assertFloatEqual(exp, obs)

        no_data = []
        exp = matrix([]).T
        obs = self.dbrda._create_factor(no_data)
        self.assertFloatEqual(exp, obs)

        mixed_data = ["foo", 20.5]
        exp = matrix([0, 1]).T
        obs = self.dbrda._create_factor(mixed_data)
        self.assertFloatEqual(exp, obs)


if __name__ == "__main__":
    main()
