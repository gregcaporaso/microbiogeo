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
from qiime.parse import parse_mapping_file_to_dict
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
    """Tests for the MetadataMap class."""

    def setUp(self):
        """Create MetadataMap objects that will be used in the tests."""
        # Create a map using the overview tutorial mapping file.
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
        self.overview_map = MetadataMap(
            *parse_mapping_file_to_dict(self.overview_map_str))

        # Create the same overview tutorial map, but this time with some
        # comments.
        self.comment = "# Some comments about this mapping file"
        self.map_with_comments_str = self.overview_map_str[:]
        self.map_with_comments_str.insert(1, self.comment)
        self.map_with_comments = MetadataMap(*parse_mapping_file_to_dict(
            self.map_with_comments_str))

        # Create a MetadataMap object that has no metadata (i.e. no sample IDs,
        # so no metadata about samples).
        self.empty_map = MetadataMap({}, [])

        # Create a MetadataMap object that has samples (i.e. sample IDs) but
        # not associated metadata (i.e. no columns other than SampleID).
        self.no_metadata_str = ["#SampleID",
                                "PC.354",
                                "PC.355",
                                "PC.356",
                                "PC.481",
                                "PC.593",
                                "PC.607",
                                "PC.634",
                                "PC.635",
                                "PC.636"]
        self.no_metadata = MetadataMap(*parse_mapping_file_to_dict(
            self.no_metadata_str))

    def test_getComments(self):
        """Test metadata map comments accessor."""
        self.assertEqual(self.overview_map.getComments(), [])
        exp = self.comment[1:]
        self.assertEqual(self.map_with_comments.getComments(), [exp])
        self.assertEqual(self.empty_map.getComments(), [])
        self.assertEqual(self.no_metadata.getComments(), [])

    def test_getSampleMetadata(self):
        """Test metadata by sample ID accessor with valid sample IDs."""
        exp = {'BarcodeSequence': 'AGCACGAGCCTA', 'Treatment': 'Control',
               'DOB': '20061218'}
        obs = self.overview_map.getSampleMetadata('PC.354')
        self.assertEqual(obs, exp)

        exp = {'BarcodeSequence': 'ACCAGCGACTAG', 'Treatment': 'Control',
               'DOB': '20070314'}
        obs = self.map_with_comments.getSampleMetadata('PC.481')
        self.assertEqual(obs, exp)

        exp = {'BarcodeSequence': 'ACGGTGAGTGTC', 'Treatment': 'Fast',
               'DOB': '20080116'}
        obs = self.map_with_comments.getSampleMetadata('PC.636')
        self.assertEqual(obs, exp)

        exp = {}
        obs = self.no_metadata.getSampleMetadata('PC.636')
        self.assertEqual(obs, exp)

    def test_getSampleMetadata_bad_sample_id(self):
        """Test metadata by sample ID accessor with invalid sample IDs."""
        # Nonexistent sample ID.
        self.assertRaises(KeyError, self.overview_map.getSampleMetadata,
            'PC.000')
        self.assertRaises(KeyError, self.no_metadata.getSampleMetadata,
            'PC.000')
        # Integer sample ID.
        self.assertRaises(KeyError, self.overview_map.getSampleMetadata, 42)
        # Sample ID of type None.
        self.assertRaises(KeyError, self.overview_map.getSampleMetadata, None)

        # Sample ID on empty map.
        self.assertRaises(KeyError, self.empty_map.getSampleMetadata, 's1')
        # Integer sample ID on empty map.
        self.assertRaises(KeyError, self.empty_map.getSampleMetadata, 1)
        # Sample ID of None on empty map.
        self.assertRaises(KeyError, self.empty_map.getSampleMetadata, None)

    def test_getCategoryValue(self):
        """Test category value by sample ID/category name accessor."""
        exp = "Fast"
        obs = self.overview_map.getCategoryValue('PC.634', 'Treatment')
        self.assertEqual(obs, exp)

        exp = "20070314"
        obs = self.overview_map.getCategoryValue('PC.481', 'DOB')
        self.assertEqual(obs, exp)

        exp = "ACGGTGAGTGTC"
        obs = self.map_with_comments.getCategoryValue(
                'PC.636', 'BarcodeSequence')
        self.assertEqual(obs, exp)

    def test_getCategoryValue_bad_sample_id(self):
        """Test category value by sample ID accessor with bad sample IDs."""
        # Nonexistent sample ID.
        self.assertRaises(KeyError, self.overview_map.getCategoryValue,
            'PC.000', 'Treatment')
        self.assertRaises(KeyError, self.no_metadata.getCategoryValue,
            'PC.000', 'Treatment')
        # Integer sample ID.
        self.assertRaises(KeyError, self.overview_map.getCategoryValue, 42,
            'DOB')
        # Sample ID of type None.
        self.assertRaises(KeyError, self.overview_map.getCategoryValue, None,
            'Treatment')

        # Sample ID on empty map.
        self.assertRaises(KeyError, self.empty_map.getCategoryValue, 's1',
            'foo')
        # Integer sample ID on empty map.
        self.assertRaises(KeyError, self.empty_map.getCategoryValue, 1,
            'bar')
        # Sample ID of None on empty map.
        self.assertRaises(KeyError, self.empty_map.getCategoryValue, None,
            'baz')

    def test_getCategoryValue_bad_category(self):
        """Test category value by sample ID accessor with bad categories."""
        # Nonexistent category.
        self.assertRaises(KeyError, self.overview_map.getCategoryValue,
            'PC.354', 'foo')
        # Integer category.
        self.assertRaises(KeyError, self.overview_map.getCategoryValue,
            'PC.354', 42)
        # Category of type None.
        self.assertRaises(KeyError, self.overview_map.getCategoryValue,
            'PC.354', None)

        # Category on map with no metadata, but that has sample IDs.
        self.assertRaises(KeyError, self.no_metadata.getCategoryValue,
            'PC.354', 'Treatment')
        # Integer category on map with no metadata.
        self.assertRaises(KeyError, self.no_metadata.getCategoryValue,
            'PC.354', 34)
        # Category of type None on map with no metadata.
        self.assertRaises(KeyError, self.no_metadata.getCategoryValue,
            'PC.354', None)


if __name__ == "__main__":
    main()
