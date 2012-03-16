#!/usr/bin/env python
# File created on 11 Mar 2012
from __future__ import division

__author__ = "Michael Dwan"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Michael Dwan", "Logan Knecht", "Jai Rideout"]
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

class GradientStatsTests(TestCase):
    """ Tests for the CategoryStats class """

    def setUp(self):
        self.test_inst = GradientStats()
  
    def test_setMetadataMap_no_instantiate(self):
        """ GradientStats is non-instantiable. Cannot call setMetadataMap method """
        self.assertRaises(NotImplementedError, self.test_inst.setMetadataMap, None)
  
    def test_getMetadataMap_no_instantiate(self):
        """ GradientStats is non-instantiable. Cannot call getMetadataMap method """
        self.assertRaises(NotImplementedError, self.test_inst.getMetadataMap)
  
    def test_setDistanceMatrix_no_instantiate(self):
        """ GradientStats is non-instantiable. Cannot call setDistanceMatrix method """
        self.assertRaises(NotImplementedError, self.test_inst.setDistanceMatrix, None)
  
    def test_getDistanceMatrix_no_instantiate(self):
        """ GradientStats is non-instantiable. Cannot call getDistanceMatrix method """
        self.assertRaises(NotImplementedError, self.test_inst.getDistanceMatrix)
  
    def test_setCategories_no_instantiate(self):
        """ GradientStats is non-instantiable. Cannot call setCategories method """
        self.assertRaises(NotImplementedError, self.test_inst.setCategories, None)

    def test_getCategories_no_instantiate(self):
        """ GradientStats is non-instantiable. Cannot call getCategories method """
        self.assertRaises(NotImplementedError, self.test_inst.getCategories)
  
    def test_runAnalysis_no_instantiate(self):
        """ GradientStats is non-instantiable. Cannot call runAnalysis method """
        self.assertRaises(NotImplementedError, self.test_inst.runAnalysis)
  

class CategoryStatsTests(TestCase):
    """ Tests for the CategoryStats class """

    def setUp(self):
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
        """ runAnalysis not implemented in abstract base CateogoryStats """ 
        raise NotImplementedError("Method no implemented by abstract base.")


class CorrelationStatsTests(TestCase):
    """ Tests for the CategoryStats class """
    #Author - LK
    def setUp(self):
        self.test_inst = CorrelationStats()

    def test_defaultVariableDeclarations(self):
        """ Basic check to make sure that the constructor is declaring the right variables with the right values """
	      #tested and verified - LK 03/13/2012
        self.assertEqual([], self.test_inst._distmat, "The _distmat default is not an empty list")

    def test_distmatGetAndSet(self):
        """ Operates on the assumption that the object was just created, tests if the getter returns the right value, an empty list """  
        #tested and verified - LK 03/13/2012
        self.test_inst = CorrelationStats()
        self.assertEqual([], self.test_inst.getDistanceMatrices(), "The default _distmat returned from getDistanceMatrices() is not an empty list")

        self.test_inst.setDistanceMatrices([1, 2, 3, 4])
        self.assertEqual([1, 2, 3, 4], self.test_inst.getDistanceMatrices(), "The setDistanceMatrices list did not set the list correctly, and the returned value is different than anticipated")

    def test_runAnalysis(self):
        """ runAnalysis not implemented in base CorrelationStats"""
	      #tested and verified - LK 03/13/2012
        self.assertRaises(NotImplementedError, self.test_inst.runAnalysis)


class DistanceMatrixStatsTests(TestCase):
    """ Tests for the DistanceMatrixStats class"""
    #Author - LK
    def setUp(self):
        self.test_inst = DistanceMatrixStats()
    
    def test_defaultVariableDeclarations(self):
        """ Basic check to make sure that the constructor is declaring the right variables with the right values """
	      #tested and verified - LK 03/13/2012
        self.assertEqual([], self.test_inst._distmat, "The _distmat default is not an empty list")

    def test_runAnalysis(self):
        """ runAnalysis not implemented in abstract base DistanceMatrixStats""" 
        self.assertRaises(NotImplementedError, self.test_inst.runAnalysis)


class MantelCorrelogramTests(TestCase):
    """Tests for the MantelCorrelogram class."""

    def setUp(self):
        # The distance matrix from the overview tutorial.
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
        self.mc.setDistanceMatrices([self.overview_dm, self.overview_dm])
        self.assertEqual(len(self.mc.getDistanceMatrices()), 2)

    def test_setDistanceMatrices_invalid(self):
        """Test setting an invalid number of distance matrices."""
        self.assertRaises(ValueError, self.mc.setDistanceMatrices,
            [self.overview_dm])
        self.assertRaises(ValueError, self.mc.setDistanceMatrices,
            [self.overview_dm, self.overview_dm, self.overview_dm])

    def test_runAnalysis(self):
        """Test running a Mantel correlogram analysis on valid input."""
        pass

if __name__ == "__main__":
    main()
