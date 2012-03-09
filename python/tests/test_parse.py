#!/usr/bin/env python

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jr378@nau.edu"
__status__ = "Development"

"""Tests functions and classes in the parse module."""
 
from numpy import array
from cogent.util.unit_test import TestCase, main
from biom.table import __version__ as __biom_version__, __url__ as __biom_url__
from python.qiime.parse import DistanceMatrix, MetadataMap

class DistanceMatrixTests(TestCase):
    """Tests for the DistanceMatrix class."""

    def setUp(self):
        """Create distance matrices that will be used in many of the tests."""
        # Create a 1x1 matrix.
        self.single_ele_dm = DistanceMatrix(array([[0]]), ['s1'], ['s1'])
        # Create a 3x3 matrix.
        self.dm = DistanceMatrix(array([[0, 2, 4], [1, 2, 3], [4, 5, 6]]),
            ['s1', 's2', 's3'], ['s1', 's2', 's3'])

    def test_getSize(self):
        """Test returning of dm's size."""
        self.assertEqual(self.single_ele_dm.getSize(), 1)
        self.assertEqual(self.dm.getSize(), 3)

    def test_empty_dm(self):
        """Can't create a dm with no data (must be at least 1x1)."""
        self.assertRaises(ValueError, DistanceMatrix, array([[]]), [], ['s1'])

    def test_nonsquare_dm(self):
        """Can't create a dm that isn't square."""
        self.assertRaises(ValueError, DistanceMatrix,
            array([[1, 2], [2, 2], [7, 4]]), ['s1', 's2'], ['s1', 's2', 's3'])

    def test_nonmatching_row_col_labels(self):
        """Can't create a dm that doesn't have matching row/col labels."""
        self.assertRaises(ValueError, DistanceMatrix, array([[1, 2], [2, 2]]),
            ['s1', 's2'], ['s1', 's3'])

    def test_nonnumpy_data(self):
        """Can't create a dm that isn't given a numpy array as data."""
        self.assertRaises(AttributeError, DistanceMatrix, [[1, 2], [2, 2]],
            ['s1', 's2'], ['s1', 's3'])

    def test_biom_type(self):
        """Make sure the BIOM type is right."""
        self.assertEqual(self.single_ele_dm._biom_type, "Distance matrix")
        self.assertEqual(self.dm._biom_type, "Distance matrix")

    def test_biom_matrix_type(self):
        """Make sure the BIOM matrix type is right."""
        self.assertEqual(self.single_ele_dm._biom_matrix_type, "dense")
        self.assertEqual(self.dm._biom_matrix_type, "dense")

    def test_getBiomFormatObject(self):
        """Should return a dictionary of the dm in BIOM format."""
        exp = {'rows': [{'id': 's1', 'metadata': None}],
               'format': 'Biological Observation Matrix %s' % __biom_version__,
               'generated_by': 'foo',
               'data': [[0]],
               'columns': [{'id': 's1', 'metadata': None}],
               'matrix_type': 'dense',
               'shape': [1, 1],
               'format_url': __biom_url__,
               'type': 'Distance matrix',
               'id': None,
               'matrix_element_type': 'int'}
        obs = self.single_ele_dm.getBiomFormatObject("foo")

        # Remove keys that we don't want to test because they might change
        # frequently (and the date is impossible to test). By using 'del', this
        # also tests that the key exists.
        del obs['date']
        self.assertEqual(obs, exp)

        exp = {'rows': [{'id': 's1', 'metadata': None},
                        {'id': 's2', 'metadata': None},
                        {'id': 's3', 'metadata': None}],
               'format': 'Biological Observation Matrix %s' % __biom_version__,
               'generated_by': 'foo',
               'data': [[0, 2, 4], [1, 2, 3], [4, 5, 6]],
               'columns': [{'id': 's1', 'metadata': None},
                           {'id': 's2', 'metadata': None},
                           {'id': 's3', 'metadata': None}],
               'matrix_type': 'dense',
               'shape': [3, 3],
               'format_url': __biom_url__,
               'type': 'Distance matrix',
               'id': None,
               'matrix_element_type': 'int'}
        obs = self.dm.getBiomFormatObject("foo")

        # Remove keys that we don't want to test because they might change
        # frequently (and the date is impossible to test). By using 'del', this
        # also tests that the key exists.
        del obs['date']
        self.assertEqual(obs, exp)


class MetadataMapTests(TestCase):
    pass


if __name__ == "__main__":
    main()
