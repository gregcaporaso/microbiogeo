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
        # The weighted unifrac distance matrix from the overview tutorial.
        self.overview_dm_str = ["\tPC.354\tPC.355\tPC.356\tPC.481\tPC.593\
                                 \tPC.607\tPC.634\tPC.635\tPC.636",
                                 "PC.354\t0.0\t0.625\t0.623\t0.61\t0.577\
                                 \t0.729\t0.8\t0.721\t0.765",
                                 "PC.355\t0.625\t0.0\t0.615\t0.642\t0.673\
                                 \t0.776\t0.744\t0.749\t0.677",
                                 "PC.356\t0.623\t0.615\t0.0\t0.682\t0.737\
                                 \t0.734\t0.777\t0.733\t0.724",
                                 "PC.481\t0.61\t0.642\t0.682\t0.0\t0.704\
                                 \t0.696\t0.675\t0.654\t0.696",
                                 "PC.593\t0.577\t0.673\t0.737\t0.704\t0.0\
                                 \t0.731\t0.758\t0.738\t0.737",
                                 "PC.607\t0.729\t0.776\t0.734\t0.696\t0.731\
                                 \t0.0\t0.718\t0.666\t0.727",
                                 "PC.634\t0.8\t0.744\t0.777\t0.675\t0.758\
                                 \t0.718\t0.0\t0.6\t0.578",
                                 "PC.635\t0.721\t0.749\t0.733\t0.654\t0.738\
                                 \t0.666\t0.6\t0.0\t0.623",
                                 "PC.636\t0.765\t0.677\t0.724\t0.696\t0.737\
                                 \t0.727\t0.578\t0.623\t0.0"]
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
        self.mc.runAnalysis()


if __name__ == "__main__":
    main()
