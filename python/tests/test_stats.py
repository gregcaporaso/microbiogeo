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
from python.qiime.stats import (BioEnv, DistanceBasedRda)

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


class NonRandomShuffler(object):
    """Helper class for testing p-values that are calculated by permutations.

    Since p-values rely on randomness, it may be useful to use a non-random
    function (such as that provided by this class) to generate permutations
    so that p-values can be accurately tested.

    This code is heavily based on Andrew Cochran's original version.
    """

    def __init__(self):
        """Default constructor initializes the number of calls to zero."""
        self.num_calls = 0

    def permutation(self, x):
        """Non-random permutation function to test p-test code.

        Returns the 'permuted' version of x.

        Arguments:
            x - the array to be 'permuted'
        """
        x = array(x)
        x = roll(x, self.num_calls)
        self.num_calls += 1
        return x


class BioEnvTests(TestHelper):
    """Tests for the BioEnv class."""

    def setUp(self):
        """Define some useful data to use in testing."""
        super(BioEnvTests, self).setUp()

        self.bv_dm_88soils_str = ["\tMT2.141698\tCA1.141704\tBB2.141659\tCO2.141657\tTL3.141709\tSN3.141650",
        "MT2.141698\t0.0\t0.623818643706\t0.750015427505\t0.585201193913\t0.729023583672\t0.622135587669",
        "CA1.141704\t0.623818643706\t0.0\t0.774881224555\t0.649822398416\t0.777203137034\t0.629507320436",
        "BB2.141659\t0.750015427505\t0.774881224555\t0.0\t0.688845424001\t0.567470311282\t0.721707516043",
        "CO2.141657\t0.585201193913\t0.649822398416\t0.688845424001\t0.0\t0.658853575764\t0.661223617505",
        "TL3.141709\t0.729023583672\t0.777203137034\t0.567470311282\t0.658853575764\t0.0\t0.711173405838",
        "SN3.141650\t0.622135587669\t0.629507320436\t0.721707516043\t0.661223617505\t0.711173405838\t0.0"]
        self.bv_dm_88soils = DistanceMatrix.parseDistanceMatrix(self.bv_dm_88soils_str)

        self.bv_map_88soils_str = ["#SampleId\tTOT_ORG_CARB\tSILT_CLAY\tELEVATION\tSOIL_MOISTURE_DEFICIT\tCARB_NITRO_RATIO\tANNUAL_SEASON_TEMP\tANNUAL_SEASON_PRECPT\tPH\tCMIN_RATE\tLONGITUDE\tLATITUDE",
        "MT2.141698\t39.1\t35\t1000\t70\t23.087\t7\t450\t6.66\t19.7\t-114\t46.8",
        "CA1.141704\t16.7\t73\t2003\t198\t13\t10.3\t400\t7.27\t2.276\t-111.7666667\t36.05",
        "BB2.141659\t52.2\t44\t400\t-680\t21.4\t6.1\t1200\t4.6\t2.223\t-68.1\t44.86666667",
        "CO2.141657\t18.1\t24\t2400\t104\t31.8\t6.1\t350\t5.68\t9.223\t-105.3333333\t40.58333333",
        "TL3.141709\t53.9\t52\t894\t-212\t24.6\t-9.3\t400\t4.23\t16.456\t-149.5833333\t68.63333333",
        "SN3.141650\t16.6\t20\t3000\t-252\t13.9\t3.6\t600\t5.74\t6.289\t-118.1666667\t36.45"]
        self.bv_map_88soils = MetadataMap.parseMetadataMap(self.bv_map_88soils_str)

        self.cats = ['TOT_ORG_CARB', 'SILT_CLAY', 'ELEVATION', 'SOIL_MOISTURE_DEFICIT', 'CARB_NITRO_RATIO', 'ANNUAL_SEASON_TEMP', 'ANNUAL_SEASON_PRECPT', 'PH', 'CMIN_RATE', 'LONGITUDE', 'LATITUDE']

        self.bioenv = BioEnv(self.bv_dm_88soils, self.bv_map_88soils, self.cats)

        self.a = [1,2,4,3,1,6,7,8,10,4]
        self.b = [2,10,20,1,3,7,5,11,6,13]
        self.c = [7,1,20,13,3,57,5,121,2,9]
        self.x = (1, 2, 4, 3, 1, 6, 7, 8, 10, 4, 100, 2, 3, 77)
        self.y = (2, 10, 20, 1, 3, 7, 5, 11, 6, 13, 5, 6, 99, 101)
        self.r = (1.7,10,20,1.7,3,7,5,11,6.5,13)
        self.s = (2,3,5,4,2,2,3,4,3,2)
        self.u = (1,2,3,4,5,6,7,8,9)
        self.v = (10,11,4,2,9,33,1,5,88)

    def test_get_rank(self):
        """Test the _get_rank method with valid input"""
        exp = ([1.5,3.5,7.5,5.5,1.5,9.0,10.0,11.0,12.0,7.5,14.0,3.5,5.5,13.0], 4)
        obs = self.bioenv._get_rank(self.x)
        self.assertFloatEqual(exp,obs)

        exp = ([1.5,3.0,5.5,4.0,1.5,7.0,8.0,9.0,10.0,5.5],2)
        obs = self.bioenv._get_rank(self.a)
        self.assertFloatEqual(exp,obs)

        exp = ([2,7,10,1,3,6,4,8,5,9],0)
        obs = self.bioenv._get_rank(self.b)
        self.assertFloatEqual(exp,obs)

        exp = ([1.5,7.0,10.0,1.5,3.0,6.0,4.0,8.0,5.0,9.0], 1)
        obs = self.bioenv._get_rank(self.r)
        self.assertFloatEqual(exp,obs)

        exp = ([],0)
        obs = self.bioenv._get_rank([])
        self.assertEqual(exp,obs)

    def test_get_rank_invalid_input(self):
        """Test the _get_rank method with invalid input"""
        vec = [1, 'a', 3, 2.5, 3, 1]
        self.assertRaises(TypeError, self.bioenv._get_rank, vec)

        vec = [1, 2, {1:2}, 2.5, 3, 1]
        self.assertRaises(TypeError, self.bioenv._get_rank, vec)

        vec = [1, 2, [23,1], 2.5, 3, 1]
        self.assertRaises(TypeError, self.bioenv._get_rank, vec)

        vec = [1, 2, (1,), 2.5, 3, 1]
        self.assertRaises(TypeError, self.bioenv._get_rank, vec)

    def test_spearman_correlation(self):
        """Test the _spearman_correlation method."""

        # One vector has no ties
        exp = 0.3719581
        obs = self.bioenv._spearman_correlation(self.a,self.b)
        self.assertFloatEqual(exp,obs)

        # Both vectors have no ties
        exp = 0.2969697
        obs = self.bioenv._spearman_correlation(self.b,self.c)
        self.assertFloatEqual(exp,obs)

        # Both vectors have ties
        exp = 0.388381
        obs = self.bioenv._spearman_correlation(self.a,self.r)
        self.assertFloatEqual(exp,obs)

    def test_spearman_correlation_invalid_input(self):
        """Test the _spearman_correlation method with invalid input."""
        self.assertRaises(ValueError,
                          self.bioenv._spearman_correlation, [],[])

        self.assertRaises(ValueError,
                          self.bioenv._spearman_correlation, self.a,[])

        self.assertRaises(ValueError,
                          self.bioenv._spearman_correlation,
                          {0:2}, [1,2,3])

    def test_vector_dist(self):
        """Test the _vector_dist helper method"""
        v1 = [1,4,2]
        v2 = [-1,12,4]

        exp = 8.48528137424
        obs = self.bioenv._vector_dist(v1,v2)
        self.assertFloatEqual(exp, obs)

        v1 = [1,2,100,4,2]
        v2 = [-1,12,4,12,99]

        exp = 137.087563258
        obs = self.bioenv._vector_dist(v1,v2)
        self.assertFloatEqual(exp, obs)

    def test_make_cat_mat(self):
        """Test the _make_cat_mat method"""
        exp = [('6.66', '46.8'),
               ('7.27', '36.05'),
               ('4.6', '44.86666667'),
               ('5.68', '40.58333333'),
               ('4.23', '68.63333333'),
               ('5.74', '36.45')]

        obs = self.bioenv._make_cat_mat(['PH', 'LONGITUDE', 'LATITUDE'], [1,3])
        self.assertEqual(exp,obs)

        exp = [('6.66', '-114', '46.8'),
               ('7.27', '-111.7666667', '36.05'),
               ('4.6', '-68.1', '44.86666667'),
               ('5.68', '-105.3333333', '40.58333333'),
               ('4.23', '-149.5833333', '68.63333333'),
               ('5.74', '-118.1666667', '36.45')]

        obs = self.bioenv._make_cat_mat(['PH', 'LONGITUDE', 'LATITUDE'], [1,2,3])
        self.assertEqual(exp,obs)

    def test_derive_euclidean_dm(self):
        """Test the derive_euclidean_dm method."""
        cat_mat = [('6.66', '46.8'),
                   ('7.27', '36.05'),
                   ('4.6', '44.86666667'),
                   ('5.68', '40.58333333'),
                   ('4.23', '68.63333333'),
                   ('5.74', '36.45')]

        exp = 0

        dm_lbls = ['MT2.141698','CA1.141704','BB2.141659',
            'CO2.141657','TL3.141709','SN3.141650']

        mtx = [
        [0.0,10.7672930674,2.82513322958,6.29343661968,
            21.9681438519,10.3908084382],
        [10.7672930674,0.0,9.21208506093,4.80408275125,
        32.7248408842,1.58142340946],
        [2.82513322958,9.21208506093,0.0,4.41739114202,
        23.7695465697,8.49351975531],
        [6.29343661968,4.80408275125,4.41739114202,0.0,
        28.0874527147,4.13376879093],
        [21.9681438519,32.7248408842,23.7695465697,
        28.0874527147,0.0,32.2187374711],
        [10.3908084382,1.58142340946,8.49351975531,
        4.13376879093,32.2187374711,0.0]]

        dm = DistanceMatrix(asarray(mtx), dm_lbls, dm_lbls)


        #print self.bioenv._derive_euclidean_dm(cat_mat,
        #                                       self.bv_dm_88soils.Size)



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
