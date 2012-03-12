#!/usr/bin/env python
# File created on 11 Mar 2012
from __future__ import division

__author__ = "Michael Dwan"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Michael Dwan, Logan Knecht"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Michael Dwan"
__email__ = "mgd25@gmail.com"
__status__ = "Development"

"""Test suite for classes, methods and functions of the stats module."""
 
from numpy import array

from cogent.util.unit_test import TestCase, main
from python.qiime.stats import GradientStats, DistanceMatrixStats, \
                               CorrelationStats, CategoryStats
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




if __name__ == "__main__":
    main()
