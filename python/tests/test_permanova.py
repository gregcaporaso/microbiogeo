#!/usr/bin/env python
# File created on 18 Jul 2011

__author__ = "Andrew Cochran"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Andrew Cochran","Dan Knights"]
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "Dan Knights"
__email__ = "danknights@gmail.com"
__status__ = "Development"

from python.qiime.permanova import permanova, _compute_f_value, permanova_p_test
from cogent.util.unit_test import TestCase, main
from python.qiime.parse import parse_distmat, parse_mapping_file_to_dict
from numpy import array, roll


class permanovaTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        self.distmtx_txt = """\tsam1\tsam2\tsam3\tsam4
        sam1\t0\t1\t5\t4
        sam2\t1\t0\t3\t2
        sam3\t5\t3\t0\t3
        sam4\t4\t2\t3\t0""".split('\n')
        self.distmtx_tie_txt = """\tsam1\tsam2\tsam3\tsam4
        sam1\t0\t1\t1\t4
        sam2\t1\t0\t3\t2
        sam3\t5\t3\t0\t3
        sam4\t4\t2\t3\t0""".split('\n')
        self.distmtx_non_sym = """\tsam1\tsam2\tsam3\tsam4\tsam5
        sam1\t0\t3\t7\t2\t1
        sam2\t3\t0\t5\t4\t1
        sam3\t7\t5\t0\t2\t6
        sam4\t2\t4\t2\t0\t2
        sam5\t1\t1\t6\t6\t0""".split('\n')
        self.mapping_txt="""#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tTreatment\tDOB\tDescription
        #Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
        sam1\tAGCACGAGCCTA\tYATGCTGCCTCCCGTAGGAGT\tControl\t20061218\tControl_mouse_I.D._354
        sam2\tAACTCGTCGATG\tYATGCTGCCTCCCGTAGGAGT\tControl\t20061218\tControl_mouse_I.D._355
        sam3\tACAGACCACTCA\tYATGCTGCCTCCCGTAGGAGT\tFast\t20061126\tControl_mouse_I.D._356
        sam4\tACCAGCGACTAG\tYATGCTGCCTCCCGTAGGAGT\tFast\t20070314\tControl_mouse_I.D._481""".split('\n')
        self.mapping_non_sym="""#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tTreatment\tDOB\tDescription
        #Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
        sam1\tAGCACGAGCCTA\tYATGCTGCCTCCCGTAGGAGT\tControl\t20061218\tControl_mouse_I.D._354
        sam2\tAACTCGTCGATG\tYATGCTGCCTCCCGTAGGAGT\tControl\t20061218\tControl_mouse_I.D._355
        sam3\tACAGACCACTCA\tYATGCTGCCTCCCGTAGGAGT\tFast\t20061126\tControl_mouse_I.D._356
        sam4\tACCAGCGACTAG\tYATGCTGCCTCCCGTAGGAGT\tAwesome\t20070314\tControl_mouse_I.D._481
        sam5\tACCAGCGACTAG\tYATGCTGCCTCCCCTATADST\tAwesome\t202020\tcontrolmouseid""".split('\n')


    def test_permanova1(self):
        """permanova should return 4.4"""
        group_list = {}
        samples, distmtx = parse_distmat(self.distmtx_txt)
        dict, comment = parse_mapping_file_to_dict(self.mapping_txt)
        for sample in dict:
            group_list[sample] = dict[sample]["Treatment"]
        result = permanova(samples, distmtx, group_list)
        self.assertEqual(result, 4.4)
        
    def test_permanova2(self):
        """Should result in 2"""
        group_list = {}
        samples, distmtx = parse_distmat(self.distmtx_tie_txt)
        dict, comment = parse_mapping_file_to_dict(self.mapping_txt)
        for sample in dict:
            group_list[sample] = dict[sample]["Treatment"]
        result = permanova(samples, distmtx, group_list)
        self.assertEqual(result, 2)
        
    def test_permanova3(self):
        """Should result in 3.58462"""
        group_list = {}
        samples, distmtx = parse_distmat(self.distmtx_non_sym)
        dict, comment = parse_mapping_file_to_dict(self.mapping_non_sym)
        for sample in dict:
            group_list[sample] = dict[sample]["Treatment"]
        result = permanova(samples, distmtx, group_list)
        self.assertEqual(round(result,5), 3.58462) 
        
    def test_compute_f1(self):
        """Should return 4.4, testing just function"""
        distances = [1,5,4,3,2,3]
        grouping = [0,-1,-1,-1,-1,1]
        distances = array(distances)
        grouping = array(grouping)
        result = _compute_f_value(distances,grouping,4,2,[2,2])
        self.assertEqual(result, 4.4)
        
    def test_p_test(self):
        """P-value should be .5 for this test"""
        group_list = {}
        samples, distmtx = parse_distmat(self.distmtx_txt)
        grouping, comment = parse_mapping_file_to_dict(self.mapping_txt)
        for sample in grouping:
            group_list[sample] = grouping[sample]["Treatment"]

        nrs = NonRandomShuffler()
	
	print(group_list)
        result, p_val = permanova_p_test(samples, distmtx,group_list, ntrials=3, randomfun=nrs.permutation)
        self.assertEqual(p_val, 0.5)


class NonRandomShuffler():
    """
    Class for testing the p_test function
        Since the p_test function relies on randomness, it is necessary to use
        a function as one of the parameters so that it can be accuratly tested
    """
    def __init__(self):
        self.numcalls = 0

    def permutation(self, x):
        """
        Non-random permutation function to test p-test code

            PARAMETERS
            self: in order to use the numcalls variable
            x: the array to be 'permutated'

            RETURNS
            x: the 'permuted' array
        """

        x = array(x)
        x = roll(x,self.numcalls)
        self.numcalls += 1        
        return x


#run tests if called from command line
if __name__ == '__main__':
    main()
