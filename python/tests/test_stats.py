#!/usr/bin/env python
# File created on 11 Mar 2012
from __future__ import division

__author__ = "Michael Dwan"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Michael Dwan", "Logan Knecht", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Michael Dwan"
__email__ = "mdwan.tgen@gmail.com"
__status__ = "Development"

"""Test suite for classes, methods and functions of the stats module."""
 
from numpy import array
from math import sqrt

from cogent.util.unit_test import TestCase, main
from python.qiime.stats import GradientStats, DistanceMatrixStats, \
                               CorrelationStats, CategoryStats, \
                               MantelCorrelogram, PartialMantel
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

        # Mantel correlogram test using the overview tutorial's unifrac dm as
        # both inputs.
        self.mc = MantelCorrelogram(self.overview_dm, self.overview_dm, 999)

        # Smallest test case: 3x3 matrices.
        ids = ['s1', 's2', 's3']
        self.small_mc = MantelCorrelogram(
            DistanceMatrix(array([[0, 1, 2], [1, 0, 3], [2, 3, 0]]), ids, ids),
            DistanceMatrix(array([[0, 2, 5], [2, 0, 8], [5, 8, 0]]), ids, ids),
            999)
    
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

    def test_getAlpha(self):
        """Test retrieving the value of alpha."""
        self.assertEqual(self.mc.getAlpha(), 0.05)

    def test_setAlpha(self):
        """Test setting the value of alpha."""
        self.mc.setAlpha(0.01)
        self.assertEqual(self.mc.getAlpha(), 0.01)

    def test_setAlpha(self):
        """Test setting the value of alpha with an invalid value."""
        self.assertRaises(ValueError, self.mc.setAlpha, -5)
        self.assertRaises(ValueError, self.mc.setAlpha, 2)

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

    def test_setDistanceMatrices_too_small(self):
        """Test setting distance matrices that are too small."""
        self.assertRaises(ValueError, self.mc.setDistanceMatrices,
            [self.single_ele_dm, self.single_ele_dm])

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
        self.assertTrue(p_vals[0] >= 0 and p_vals[0] <= 0.01)
        self.assertTrue(p_vals[1] > 0.01 and p_vals[1] <= 0.1)
        self.assertTrue(p_vals[2] > 0.1 and p_vals[2] <= 0.5)
        self.assertEqual(p_vals[3:], [None, None, None, None])
        self.assertFloatEqual(corr_p_vals,
            [p_val * 3 if p_val is not None else None for p_val in p_vals])

    def test_runAnalysis_small(self):
        """Test running a Mantel correlogram analysis on the smallest input."""
        # The expected output was verified with vegan's mantel correlogram
        # function.
        obs = self.small_mc.runAnalysis()

        exp_method_name = 'Mantel Correlogram'
        self.assertEqual(obs['method_name'], exp_method_name)

        exp_class_index = [3.0, 5.0, 7.0]
        self.assertFloatEqual(obs['class_index'], exp_class_index)

        exp_num_dist = [2, 2, 2]
        self.assertEqual(obs['num_dist'], exp_num_dist)

        exp_mantel_r = [0.86602540378443871, None, None]
        self.assertFloatEqual(obs['mantel_r'], exp_mantel_r)

        # Test matplotlib Figure for a sane state.
        obs_fig = obs['correlogram_plot']
        obs_ax = obs_fig.get_axes()[0]
        self.assertEqual(obs_ax.get_title(), "Mantel Correlogram")
        self.assertEqual(obs_ax.get_xlabel(), "Distance class index")
        self.assertEqual(obs_ax.get_ylabel(), "Mantel correlation statistic")
        self.assertFloatEqual(obs_ax.get_xticks(), [2.85, 2.9, 2.95, 3., 3.05,
            3.1, 3.15, 3.2])
        self.assertFloatEqual(obs_ax.get_yticks(), [0.82, 0.83, 0.84, 0.85,
            0.86, 0.87, 0.88, 0.89, 0.9, 0.91])

        # Test p-values and corrected p-values.
        p_vals = obs['mantel_p']
        corr_p_vals = obs['mantel_p_corr']
        self.assertEqual(len(p_vals), 3)
        self.assertTrue(p_vals[0] >= 0 and p_vals[0] <= 0.5)
        self.assertEqual(p_vals[1:], [None, None])
        self.assertFloatEqual(corr_p_vals, p_vals)

    def test_find_distance_classes(self):
        """Test finding the distance classes a matrix's elements are in."""
        exp = (array([[-1,  0,  1], [ 0, -1,  2], [ 1,  2, -1]]),
               [3.0, 5.0, 7.0])
        obs = self.small_mc._find_distance_classes(
            self.small_mc.getDistanceMatrices()[1], 3)
        self.assertFloatEqual(obs, exp)

        exp = (array([[-1,  1,  2,  0,  0,  5,  7,  4,  6],
            [ 1, -1,  0,  2,  3,  6,  6,  6,  4],
            [ 2,  0, -1,  4,  5,  5,  7,  4,  6],
            [ 0,  2,  4, -1,  3,  3,  3,  3,  2],
            [ 0,  3,  5,  3, -1,  5,  7,  6,  6],
            [ 5,  6,  5,  3,  5, -1,  5,  2,  5],
            [ 7,  6,  7,  3,  7,  5, -1,  0,  0],
            [ 4,  6,  4,  3,  6,  2,  0, -1,  0],
            [ 6,  4,  6,  2,  6,  5,  0,  0, -1]]),
            [0.57381779, 0.60024231, 0.62666684, 0.65309137, 0.67951589,
             0.70594042, 0.73236494, 0.75878947])
        obs = self.mc._find_distance_classes(
            self.mc.getDistanceMatrices()[1], 8)
        self.assertFloatEqual(obs, exp)

    def test_find_distance_classes_invalid_num_classes(self):
        """Test finding the distance classes for a bad number of classes."""
        self.assertRaises(ValueError, self.mc._find_distance_classes,
                self.mc.getDistanceMatrices()[1], 0)
        self.assertRaises(ValueError, self.mc._find_distance_classes,
                self.mc.getDistanceMatrices()[1], -1)

    def test_find_break_points(self):
        """Test finding equal-spaced breakpoints in a range."""
        exp = [-2.2204460492503131e-16, 1.0, 2.0, 3.0, 4.0, 5.0]
        obs = self.mc._find_break_points(0, 5, 5)
        self.assertFloatEqual(obs, exp)

        exp = [-2.0, -1.66666666667, -1.33333333333, -1.0]
        obs = self.mc._find_break_points(-2, -1, 3)
        self.assertFloatEqual(obs, exp)

        exp = [-1.0, -0.5, 0.0, 0.5, 1.0]
        obs = self.mc._find_break_points(-1, 1, 4)
        self.assertFloatEqual(obs, exp)

        exp = [-1.0, 1.0]
        obs = self.mc._find_break_points(-1, 1, 1)
        self.assertFloatEqual(obs, exp)

    def test_find_break_points_invalid_range(self):
        """Test finding breakpoints on an invalid range."""
        self.assertRaises(ValueError, self.mc._find_break_points, 1, 0, 5)
        self.assertRaises(ValueError, self.mc._find_break_points, 1, 1, 5)

    def test_find_break_points_invalid_num_classes(self):
        """Test finding breakpoints with an invalid number of classes."""
        self.assertRaises(ValueError, self.mc._find_break_points, 0, 1, 0)
        self.assertRaises(ValueError, self.mc._find_break_points, 0, 1, -1)

    def test_generate_correlogram(self):
        """Test creating a correlogram plot."""
        obs_fig = self.mc._generate_correlogram([0, 1, 2], [-0.9, 0, 0.9],
                [0.001, 0.1, 0.9])
        obs_ax = obs_fig.get_axes()[0]
        self.assertEqual(obs_ax.get_title(), "Mantel Correlogram")
        self.assertEqual(obs_ax.get_xlabel(), "Distance class index")
        self.assertEqual(obs_ax.get_ylabel(), "Mantel correlation statistic")
        self.assertFloatEqual(obs_ax.get_xticks(), [0., 0.5, 1., 1.5, 2.])
        self.assertFloatEqual(obs_ax.get_yticks(), [-1., -0.5, 0., 0.5, 1.])

    def test_generate_correlogram_empty(self):
        """Test creating a correlogram plot with no data."""
        obs_fig = self.mc._generate_correlogram([], [], [])
        obs_ax = obs_fig.get_axes()[0]
        self.assertEqual(obs_ax.get_title(), "Mantel Correlogram")
        self.assertEqual(obs_ax.get_xlabel(), "Distance class index")
        self.assertEqual(obs_ax.get_ylabel(), "Mantel correlation statistic")
        self.assertFloatEqual(obs_ax.get_xticks(),
            [0., 0.2, 0.4, 0.6, 0.8, 1.0])
        self.assertFloatEqual(obs_ax.get_yticks(),
            [0., 0.2, 0.4, 0.6, 0.8, 1.0])

class MantelTests(TestHelper):
    """Tests for the Mantel class."""
    def setUp(self):
        """Set up PartialMantel instances for use in tests."""
        super(MantelTests, self).setUp()

        self.mantel = Mantel(self.overview_dm, self.overview_dm, 999)

        smpl_ids = ['s1', 's2', 's3']
        self.small_mantel = Mantel(
            DistanceMatrix(array([[1, 3, 2], [1, 1, 3], [4, 3, 1]]), smpl_ids, smpl_ids),
            DistanceMatrix(array([[0, 2, 5], [2, 0, 8], [5, 8, 0]]), smpl_ids, smpl_ids),
            999)

class PartialMantelTests(TestHelper):
    """Tests for the PartialMantel class."""

    def setUp(self):
        """Set up PartialMantel instances for use in tests."""
        super(PartialMantelTests, self).setUp()

        # Test partial Mantel using the unifrac dm from the overview tutorial as
        # all three inputs(should be a small value).
        self.pm = PartialMantel(self.overview_dm, self.overview_dm, self.overview_dm, 999)

        # Justa small matrix that is easy to edit and observe.
        smpl_ids = ['s1', 's2', 's3']
        self.small_pm = PartialMantel(
            DistanceMatrix(array([[1, 3, 2], [1, 1, 3], [4, 3, 1]]), smpl_ids, smpl_ids),
            DistanceMatrix(array([[0, 2, 5], [2, 0, 8], [5, 8, 0]]), smpl_ids, smpl_ids),
            DistanceMatrix(array([[10, 7, 13], [9, 7, 0], [10, 2, 8]]), smpl_ids, smpl_ids),
            999)
    
    def test_getNumPermutations(self):
        """Test retrieval of the number of permutations."""
        self.assertEqual(self.pm.getNumPermutations(), 999)

    def test_setNumPermutations(self):
        """Test setting of the number of permutations."""
        self.pm.setNumPermutations(7)
        self.assertEqual(self.pm.getNumPermutations(), 7)

    def test_setNumPermutations_invalid(self):
        """Test setting of the number of permutations using a negative(invalid) number."""
        self.assertRaises(ValueError, self.pm.setNumPermutations, -22)

    def test_setDistanceMatrices(self):
        """Test setting matrices using a valid number of distance matrices."""
        dms = [self.overview_dm, self.overview_dm, self.overview_dm]
        self.pm.setDistanceMatrices(dms)
        self.assertEqual(self.pm.getDistanceMatrices(), dms)

    def test_setDistanceMatrices_wrong_number(self):
        """Test setting matrices using an invalid number of distance matrices."""
        self.assertRaises(ValueError, self.pm.setDistanceMatrices,
            [self.overview_dm, self.overview_dm])
        self.assertRaises(ValueError, self.pm.setDistanceMatrices,
            [self.overview_dm, self.overview_dm, self.overview_dm, self.overview_dm])

    def test_runAnalysis(self):
        """Test running partial Mantel analysis on valid input."""
        obs = self.pm.runAnalysis()

        #print '\n' + str(obs) + '\n'

    def test_runAnalysis(self):
        """Test the correct running of partial Mantel analysis on small, controlled input."""
        # The output needs to be verified against the Vegan mantel.partial
        # function.
        obs = self.small_pm.runAnalysis()

        #print '\n' + str(obs) + '\n'


if __name__ == "__main__":
    main()
