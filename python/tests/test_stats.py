#!/usr/bin/env python
# File created on 11 Mar 2012
from __future__ import division

__author__ = "Michael Dwan"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Michael Dwan", "Logan Knecht", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Michael Dwan"
__email__ = "mgd25@gmail.com"
__status__ = "Development"

"""Test suite for classes, methods and functions of the stats module."""
 
from numpy import array

from cogent.util.unit_test import TestCase, main
from python.qiime.stats import GradientStats, DistanceMatrixStats, \
                               CorrelationStats, CategoryStats, \
                               MantelCorrelogram
from python.qiime.parse import DistanceMatrix, MetadataMap

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

        # A 1x1 dm.
        self.single_ele_dm = DistanceMatrix(array([[0]]), ['s1'], ['s1'])

class GradientStatsTests(TestCase):
    """Tests for the GradientStats class."""

    def setUp(self):
        self.test_inst = GradientStats()
  
    def test_runAnalysis_no_instantiate(self):
        """GradientStats is non-instantiable, cannot call runAnalysis()."""
        self.assertRaises(NotImplementedError, self.test_inst.runAnalysis)


class DistanceMatrixStatsTests(TestHelper):
    """Tests for the DistanceMatrixStats class."""

    def setUp(self):
        """Define some dm stats instances that will be used by the tests."""
        super(DistanceMatrixStatsTests, self).setUp()

        self.empty_dms = DistanceMatrixStats([])
        self.single_dms = DistanceMatrixStats([self.overview_dm])
        self.double_dms = DistanceMatrixStats(
            [self.overview_dm, self.single_ele_dm])
    
    def test_getDistanceMatrices(self):
        """Test getter for distmats."""
        self.assertEqual(self.empty_dms.getDistanceMatrices(), [])
        self.assertEqual(self.single_dms.getDistanceMatrices(),
            [self.overview_dm])
        self.assertEqual(self.double_dms.getDistanceMatrices(),
            [self.overview_dm, self.single_ele_dm])

    def test_setDistanceMatrices(self):
        """Test setter for dms on valid input data."""
        self.empty_dms.setDistanceMatrices([])
        self.assertEqual(self.empty_dms.getDistanceMatrices(), [])

        self.empty_dms.setDistanceMatrices([self.overview_dm])
        self.assertEqual(self.empty_dms.getDistanceMatrices(),
            [self.overview_dm])

        self.empty_dms.setDistanceMatrices(
            [self.overview_dm, self.overview_dm])
        self.assertEqual(self.empty_dms.getDistanceMatrices(),
            [self.overview_dm, self.overview_dm])

    def test_setDistanceMatrices_invalid(self):
        """Test setter for dms on invalid input data."""
        self.assertRaises(TypeError, self.empty_dms.setDistanceMatrices, None)
        self.assertRaises(TypeError, self.empty_dms.setDistanceMatrices, 10)
        self.assertRaises(TypeError, self.empty_dms.setDistanceMatrices, 20.0)
        self.assertRaises(TypeError, self.empty_dms.setDistanceMatrices, "foo")
        self.assertRaises(TypeError, self.empty_dms.setDistanceMatrices, {})
        self.assertRaises(TypeError, self.empty_dms.setDistanceMatrices,
            self.overview_dm)

        # Test constructor as well.
        self.assertRaises(TypeError, DistanceMatrixStats, None)
        self.assertRaises(TypeError, DistanceMatrixStats, 10)
        self.assertRaises(TypeError, DistanceMatrixStats, 20.0)
        self.assertRaises(TypeError, DistanceMatrixStats, "foo")
        self.assertRaises(TypeError, DistanceMatrixStats, {})
        self.assertRaises(TypeError, DistanceMatrixStats, self.overview_dm)


class CorrelationStatsTests(TestHelper):
    """Tests for the CorrelationStats class."""

    def setUp(self):
        """Set up correlation stats instances for use in tests."""
        super(CorrelationStatsTests, self).setUp()
        self.cs = CorrelationStats([self.overview_dm, self.overview_dm])

    def test_setDistanceMatrices(self):
        """Test setting valid distance matrices."""
        dms = [self.overview_dm, self.overview_dm]
        self.cs.setDistanceMatrices(dms)
        self.assertEqual(self.cs.getDistanceMatrices(), dms)

        dms = [self.overview_dm, self.overview_dm, self.overview_dm]
        self.cs.setDistanceMatrices(dms)
        self.assertEqual(self.cs.getDistanceMatrices(), dms)

    def test_setDistanceMatrices_too_few(self):
        """Test setting dms with not enough of them."""
        self.assertRaises(ValueError, self.cs.setDistanceMatrices, [])
        # Also test that constructor raises this error.
        self.assertRaises(ValueError, CorrelationStats, [])

    def test_setDistanceMatrices_wrong_dims(self):
        """Test setting dms with mismatching dimensions."""
        self.assertRaises(ValueError, self.cs.setDistanceMatrices,
            [self.overview_dm, self.single_ele_dm])
        # Also test that constructor raises this error.
        self.assertRaises(ValueError, CorrelationStats, [self.overview_dm,
                          self.single_ele_dm])

    def test_setDistanceMatrices_mismatched_labels(self):
        """Test setting dms with mismatching sample ID labels."""
        mismatch = DistanceMatrix(array([[0]]), ['s2'], ['s2'])
        self.assertRaises(ValueError, self.cs.setDistanceMatrices,
            [self.single_ele_dm, mismatch])
        # Also test that constructor raises this error.
        self.assertRaises(ValueError, CorrelationStats, [self.single_ele_dm,
                          mismatch])

    def test_runAnalysis(self):
        """Test runAnalysis() is not implemented in CorrelationStats"""
        self.assertRaises(NotImplementedError, self.cs.runAnalysis)


class CategoryStatsTests(TestCase):
    """Tests for the CategoryStats class."""

    def setUp(self):
        """Define some useful data to use in testing."""
        self.test_map = MetadataMap({}, [])
        self.test_dm = DistanceMatrix(array([[0]]), ['s1'], ['s1'])
        self.test_cats = ['cat1', 'cat2', 'cat3']
        self.test_inst = CategoryStats(self.test_map, self.test_dm, self.test_cats)

    def test_setMetadataMap_valid_input(self):
        """ CategoryStats.setMetadataMap() must receive an instance of MetadataMap """
        self.assertRaises(TypeError, self.test_inst.setMetadataMap, "Hello!")
        self.assertRaises(TypeError, self.test_inst.setMetadataMap, self.test_dm)
        #etc...

    def test_getMetadataMap(self):
        """ Test valid return of getMetadataMap method """
        expected = MetadataMap({}, [])
        observed = self.test_inst.getMetadataMap()
        self.assertEqual(expected, observed)

    def test_setDistanceMatrix_valid_input(self):
        """ CategoryStats.setDistanceMatrix() must receive an instance of DistanceMatrix """
        self.assertRaises(TypeError, self.test_inst.setDistanceMatrix, "Hello!")
        self.assertRaises(TypeError, self.test_inst.setDistanceMatrix, self.test_map)

    def test_getDistanceMatrix(self):
        """ Test valid return of getDistanceMatrix method """
        expected = DistanceMatrix(array([[0]]), ['s1'], ['s1'])
        observed = self.test_inst.getDistanceMatrix()
        self.assertEqual(expected, observed)

    def test_setCategories_valid_input(self):
        """ CategoryStats.setCategories() must receive an a list of strings """
        self.assertRaises(TypeError, self.test_inst.setCategories, "Hello!")
        self.assertRaises(TypeError, self.test_inst.setCategories, self.test_dm)
        self.assertRaises(TypeError, self.test_inst.setCategories, ["hehehe", 123, "hello"])

    def test_getCategories(self):
        """ Test valid return of getDistanceMatrix method """
        expected = ['cat1', 'cat2', 'cat3']
        observed = self.test_inst.getCategories()
        self.assertEqual(expected, observed)

    def runAnalysis(self):
        """ runAnalysis not implemented in abstract base CategoryStats """ 
        raise NotImplementedError("Method not implemented by abstract base.")


class MantelCorrelogramTests(TestHelper):
    """Tests for the MantelCorrelogram class."""

    def setUp(self):
        """Set up mantel correlogram instances for use in tests."""
        super(MantelCorrelogramTests, self).setUp()
        self.mc = MantelCorrelogram(self.overview_dm, self.overview_dm, 999)
    
    def test_getNumPermutations(self):
        """Test retrieving the number of permutations."""
        self.assertEqual(self.mc.getNumPermutations(), 999)

    def test_setNumPermutations(self):
        """Test setting the number of permutations."""
        self.mc.setNumPermutations(5)
        self.assertEqual(self.mc.getNumPermutations(), 5)

    def test_setNumPermutations_invalid(self):
        """Test setting the number of permutations with a negative number."""
        self.assertRaises(ValueError, self.mc.setNumPermutations, -5)

    def test_setDistanceMatrices(self):
        """Test setting a valid number of distance matrices."""
        dms = [self.overview_dm, self.overview_dm]
        self.mc.setDistanceMatrices(dms)
        self.assertEqual(self.mc.getDistanceMatrices(), dms)

    def test_setDistanceMatrices_wrong_number(self):
        """Test setting an invalid number of distance matrices."""
        self.assertRaises(ValueError, self.mc.setDistanceMatrices,
            [self.overview_dm])
        self.assertRaises(ValueError, self.mc.setDistanceMatrices,
            [self.overview_dm, self.overview_dm, self.overview_dm])

    def test_runAnalysis(self):
        """Test running a Mantel correlogram analysis on valid input."""
        # A lot of the returned numbers are based on random permutations and
        # thus cannot be tested for exact values. We'll test what we can
        # exactly, and then test for "sane" values for the "random" values. The
        # matplotlib Figure object cannot be easily tested either, so we'll try
        # our best to make sure it appears sane.
        obs = self.mc.runAnalysis()

        exp_method_name = 'Mantel Correlogram'
        self.assertEqual(obs['method_name'], exp_method_name)

        exp_class_index = [0.5757052546507142, 0.60590471266814283,
            0.63610417068557146, 0.66630362870299997, 0.69650308672042849,
            0.72670254473785723, 0.75690200275528574]
        self.assertFloatEqual(obs['class_index'], exp_class_index)

        exp_num_dist = [12, 6, 8, 10, 12, 16, 8]
        self.assertEqual(obs['num_dist'], exp_num_dist)

        exp_mantel_r = [0.73244729118260765, 0.31157641757444593,
            0.17627427296718071, None, None, None, None]
        self.assertFloatEqual(obs['mantel_r'], exp_mantel_r)

        # Test matplotlib Figure for a sane state.
        obs_fig = obs['correlogram_plot']
        obs_ax = obs_fig.get_axes()[0]
        self.assertEqual(obs_ax.get_title(), "Mantel Correlogram")
        self.assertEqual(obs_ax.get_xlabel(), "Distance class index")
        self.assertEqual(obs_ax.get_ylabel(), "Mantel correlation statistic")
        self.assertFloatEqual(obs_ax.get_xticks(), [0.57, 0.58, 0.59, 0.6,
            0.61, 0.62, 0.63, 0.64, 0.65])
        self.assertFloatEqual(obs_ax.get_yticks(), [0.1, 0.2, 0.3, 0.4, 0.5,
            0.6, 0.7, 0.8, 0.9])

        # Test p-values and corrected p-values.
        p_vals = obs['mantel_p']
        corr_p_vals = obs['mantel_p_corr']
        self.assertEqual(len(p_vals), 7)
        self.assertTrue(p_vals[0] <= 0.01)
        self.assertTrue(p_vals[1] > 0.01 and p_vals[1] <= 0.1)
        self.assertTrue(p_vals[2] > 0.1 and p_vals[2] <= 0.5)
        self.assertEqual(p_vals[3:], [None, None, None, None])
        self.assertFloatEqual(corr_p_vals,
            [p_val * 3 if p_val is not None else None for p_val in p_vals])
        

if __name__ == "__main__":
    main()
