#!/usr/bin/env python

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jai Ram Rideout, Logan Knecht"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jr378@nau.edu"
__status__ = "Development"

"""Tests functions and classes in the parse module."""
 
from numpy import array
from cogent.util.unit_test import TestCase, main
from biom.table import __version__ as __biom_version__, __url__ as __biom_url__
from qiime.parse import parse_distmat, parse_mapping_file_to_dict, \
    QiimeParseError
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

        # A distance matrix similar to the overview tutorial's unifrac dm. I
        # found this in some other tests in QIIME, but these values don't match
        # the values found in the overview tutorial's unweighted or weighted
        # unifrac distance matrices, so I'm not quite sure where this data came
        # from. That's okay, though, as we can still use it for the unit tests.
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
        sample_ids, matrix_data = parse_distmat(self.overview_dm_str)
        self.overview_dm = DistanceMatrix(matrix_data, sample_ids, sample_ids)

    def test_parseDistanceMatrix(self):
        """Test parsing a distance matrix into a DistanceMatrix instance."""
        obs = DistanceMatrix.parseDistanceMatrix(self.overview_dm_str)
        self.assertFloatEqual(obs, self.overview_dm)

    def test_parseDistanceMatrix_empty(self):
        """Test parsing empty dm file contents."""
        self.assertRaises(TypeError, DistanceMatrix.parseDistanceMatrix, [])

    def test_getSize(self):
        """Test returning of dm's size."""
        self.assertEqual(self.single_ele_dm.getSize(), 1)
        self.assertEqual(self.dm.getSize(), 3)
        self.assertEqual(self.overview_dm.getSize(), 9)

    def test_flatten(self):
        """Test flattening various dms."""
        self.assertEqual(self.single_ele_dm.flatten(), [])
        self.assertEqual(self.dm.flatten(), [1, 4, 5])
        exp = [0.625, 0.623, 0.60999999999999999, 0.57699999999999996,
            0.72899999999999998, 0.80000000000000004, 0.72099999999999997,
            0.76500000000000001, 0.61499999999999999, 0.64200000000000002,
            0.67300000000000004, 0.77600000000000002, 0.74399999999999999,
            0.749, 0.67700000000000005, 0.68200000000000005,
            0.73699999999999999, 0.73399999999999999, 0.77700000000000002,
            0.73299999999999998, 0.72399999999999998, 0.70399999999999996,
            0.69599999999999995, 0.67500000000000004, 0.65400000000000003,
            0.69599999999999995, 0.73099999999999998, 0.75800000000000001,
            0.73799999999999999, 0.73699999999999999, 0.71799999999999997,
            0.66600000000000004, 0.72699999999999998, 0.59999999999999998,
            0.57799999999999996, 0.623]
        self.assertFloatEqual(self.overview_dm.flatten(), exp)

    def test_flatten_all(self):
        """Test flattening various dms including all elements."""
        self.assertEqual(self.single_ele_dm.flatten(lower=False), [0])
        self.assertEqual(self.dm.flatten(lower=False),
            [0, 1, 4, 2, 2, 5, 4, 3, 6])
        exp = [0.0, 0.625, 0.623, 0.60999999999999999, 0.57699999999999996,
            0.72899999999999998, 0.80000000000000004, 0.72099999999999997,
            0.76500000000000001, 0.625, 0.0, 0.61499999999999999,
            0.64200000000000002, 0.67300000000000004, 0.77600000000000002,
            0.74399999999999999, 0.749, 0.67700000000000005, 0.623,
            0.61499999999999999, 0.0, 0.68200000000000005, 0.73699999999999999,
            0.73399999999999999, 0.77700000000000002, 0.73299999999999998,
            0.72399999999999998, 0.60999999999999999, 0.64200000000000002,
            0.68200000000000005, 0.0, 0.70399999999999996, 0.69599999999999995,
            0.67500000000000004, 0.65400000000000003, 0.69599999999999995,
            0.57699999999999996, 0.67300000000000004, 0.73699999999999999,
            0.70399999999999996, 0.0, 0.73099999999999998, 0.75800000000000001,
            0.73799999999999999, 0.73699999999999999, 0.72899999999999998,
            0.77600000000000002, 0.73399999999999999, 0.69599999999999995,
            0.73099999999999998, 0.0, 0.71799999999999997, 0.66600000000000004,
            0.72699999999999998, 0.80000000000000004, 0.74399999999999999,
            0.77700000000000002, 0.67500000000000004, 0.75800000000000001,
            0.71799999999999997, 0.0, 0.59999999999999998, 0.57799999999999996,
            0.72099999999999997, 0.749, 0.73299999999999998,
            0.65400000000000003, 0.73799999999999999, 0.66600000000000004,
            0.59999999999999998, 0.0, 0.623, 0.76500000000000001,
            0.67700000000000005, 0.72399999999999998, 0.69599999999999995,
            0.73699999999999999, 0.72699999999999998, 0.57799999999999996,
            0.623, 0.0]
        self.assertFloatEqual(self.overview_dm.flatten(lower=False), exp)

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
        self.assertEqual(self.overview_dm._biom_type, "Distance matrix")

    def test_biom_matrix_type(self):
        """Make sure the BIOM matrix type is right."""
        self.assertEqual(self.single_ele_dm._biom_matrix_type, "dense")
        self.assertEqual(self.dm._biom_matrix_type, "dense")
        self.assertEqual(self.overview_dm._biom_matrix_type, "dense")

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
        del obs['date']
        self.assertEqual(obs, exp)

        exp = {'rows': [{'id': 'PC.354', 'metadata': None}, {'id': 'PC.355',
            'metadata': None}, {'id': 'PC.356', 'metadata': None}, {'id':
            'PC.481', 'metadata': None}, {'id': 'PC.593', 'metadata': None},
            {'id': 'PC.607', 'metadata': None}, {'id': 'PC.634', 'metadata':
            None}, {'id': 'PC.635', 'metadata': None}, {'id': 'PC.636',
            'metadata': None}], 'format':
            'Biological Observation Matrix 0.9.1-dev', 'data': [[0.0, 0.625,
            0.623, 0.60999999999999999, 0.57699999999999996,
            0.72899999999999998, 0.80000000000000004, 0.72099999999999997,
            0.76500000000000001], [0.625, 0.0, 0.61499999999999999,
            0.64200000000000002, 0.67300000000000004, 0.77600000000000002,
            0.74399999999999999, 0.749, 0.67700000000000005], [0.623,
            0.61499999999999999, 0.0, 0.68200000000000005, 0.73699999999999999,
            0.73399999999999999, 0.77700000000000002, 0.73299999999999998,
            0.72399999999999998], [0.60999999999999999, 0.64200000000000002,
            0.68200000000000005, 0.0, 0.70399999999999996, 0.69599999999999995,
            0.67500000000000004, 0.65400000000000003, 0.69599999999999995],
            [0.57699999999999996, 0.67300000000000004, 0.73699999999999999,
            0.70399999999999996, 0.0, 0.73099999999999998, 0.75800000000000001,
            0.73799999999999999, 0.73699999999999999], [0.72899999999999998,
            0.77600000000000002, 0.73399999999999999, 0.69599999999999995,
            0.73099999999999998, 0.0, 0.71799999999999997, 0.66600000000000004,
            0.72699999999999998], [0.80000000000000004, 0.74399999999999999,
            0.77700000000000002, 0.67500000000000004, 0.75800000000000001,
            0.71799999999999997, 0.0, 0.59999999999999998,
            0.57799999999999996], [0.72099999999999997, 0.749,
            0.73299999999999998, 0.65400000000000003, 0.73799999999999999,
            0.66600000000000004, 0.59999999999999998, 0.0, 0.623],
            [0.76500000000000001, 0.67700000000000005, 0.72399999999999998,
            0.69599999999999995, 0.73699999999999999, 0.72699999999999998,
            0.57799999999999996, 0.623, 0.0]], 'columns': [{'id': 'PC.354',
            'metadata': None}, {'id': 'PC.355', 'metadata': None}, {'id':
            'PC.356', 'metadata': None}, {'id': 'PC.481', 'metadata': None},
            {'id': 'PC.593', 'metadata': None}, {'id': 'PC.607', 'metadata':
            None}, {'id': 'PC.634', 'metadata': None}, {'id': 'PC.635',
            'metadata': None}, {'id': 'PC.636', 'metadata': None}],
            'generated_by': 'foo', 'matrix_type': 'dense', 'shape': [9, 9],
            'format_url': 'http://biom-format.org', 'type': 'Distance matrix',
            'id': None, 'matrix_element_type': 'float'}

        obs = self.overview_dm.getBiomFormatObject("foo")
        del obs['date']
        self.assertFloatEqual(obs, exp)


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

    def test_parseMetadataMap(self):
        """Test parsing a mapping file into a MetadataMap instance."""
        obs = MetadataMap.parseMetadataMap(self.overview_map_str)
        self.assertEqual(obs, self.overview_map)

    def test_parseMetadataMap_empty(self):
        """Test parsing empty mapping file contents."""
        self.assertRaises(QiimeParseError, MetadataMap.parseMetadataMap, [])

    def test_eq(self):
        """Test whether two MetadataMap's are equal."""
        self.assertTrue(self.empty_map == MetadataMap({}, []))
        self.assertTrue(self.overview_map == MetadataMap(
            self.overview_map._metadata, self.overview_map._comments))

    def test_ne(self):
        """Test whether two MetadataMap's are not equal."""
        self.assertTrue(self.empty_map != MetadataMap({}, ["foo"]))
        self.assertTrue(self.overview_map != MetadataMap(
            self.overview_map._metadata, ["foo"]))
        self.assertTrue(self.overview_map != MetadataMap({},
            self.overview_map._comments))
        self.assertTrue(self.overview_map != self.empty_map)
        self.assertTrue(self.overview_map != self.map_with_comments)
        self.assertTrue(self.overview_map != self.no_metadata)

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
