#!/usr/bin/env python


__author__ = "Andrew Cochran"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Andrew Cochran","Dan Knights"]
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "Dan Knights"
__email__ = "danknights@gmail.com"
__status__ = "Development"

from Bio.qiime.anosim import anosim, _remove_ties, _compute_r_value, anosim_p_test
from cogent.util.unit_test import TestCase, main
from qiime.parse import parse_distmat, parse_mapping_file_to_dict
from numpy import array, roll


class AnosimTests(TestCase):
    """Tests of all anosim functions"""

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
        self.mapping_txt="""#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tTreatment\tDOB\tDescription
        #Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
        sam1\tAGCACGAGCCTA\tYATGCTGCCTCCCGTAGGAGT\tControl\t20061218\tControl_mouse_I.D._354
        sam2\tAACTCGTCGATG\tYATGCTGCCTCCCGTAGGAGT\tControl\t20061218\tControl_mouse_I.D._355
        sam3\tACAGACCACTCA\tYATGCTGCCTCCCGTAGGAGT\tFast\t20061126\tControl_mouse_I.D._356
        sam4\tACCAGCGACTAG\tYATGCTGCCTCCCGTAGGAGT\tFast\t20070314\tControl_mouse_I.D._481""".split('\n')

   
    def test_anosim(self):
        """Anosim should return .625"""
        group_list = {}
        samples, distmtx = parse_distmat(self.distmtx_txt)
        grouping, comment = parse_mapping_file_to_dict(self.mapping_txt)
        for sample in grouping:
            group_list[sample] = grouping[sample]["Treatment"]
        result = anosim(samples, distmtx,group_list)
        self.assertEqual(result, 0.625) # make sure it equals .625
   
    def test_anosim_begining_tie(self):
        """Should result in [1.5 1.5 3 4.5 4.5 6] as adjusted rank"""
        group_list = {}
        samples, distmtx = parse_distmat(self.distmtx_tie_txt)
        grouping, comment = parse_mapping_file_to_dict(self.mapping_txt)
        for sample in grouping:
            group_list[sample] = grouping[sample]["Treatment"]
        result = anosim(samples, distmtx,group_list)
        self.assertEqual(result, 0.25) # make sure it equals .25
   
    def test_remove_ties1(self):
        """Should return [1.5,1.5]"""
        result = _remove_ties([1,1],[1,2])
        self.assertEqual(result, [1.5,1.5])
   
    def test_remove_ties2(self):
        """Should return [3.5,3.5,3.5,3.5,3.5,3.5]"""
        result = _remove_ties([1,1,1,1,1,1],[1,2,3,4,5,6])
        self.assertEqual(result, [3.5,3.5,3.5,3.5,3.5,3.5])
   
    def test_remove_ties3(self):
        """Should return [1,3.5,3.5,3.5,3.5,6]"""
        result = _remove_ties([1,3,3,3,3,8],[1,2,3,4,5,6])
        self.assertEqual(result, [1,3.5,3.5,3.5,3.5,6])
   
    def test_remove_ties4(self):
        """Should return [1,2,3,4]"""    
        result = _remove_ties([1,2,3,4],[1,2,3,4])
        self.assertEqual(result,[1,2,3,4])
   
    def test_remove_ties5(self):
        """Should return [1,3,3,3,5.5,5.5,7]"""
        result = _remove_ties([1,2,2,2,3,3,5],[1,2,3,4,5,6,7])
        self.assertEqual(result,[1,3,3,3,5.5,5.5,7])
   
    def test_remove_ties6(self):
        """Should return [1.5,1.5,3.5,3.5]"""
        result = _remove_ties([1,1,2,2],[1,2,3,4])
        self.assertEqual(result,[1.5,1.5,3.5,3.5])
   
    def test_compute_r1(self):
        """Should return .625"""
        sorted_rank = [1.0,2.0,3.5,3.5,5.0,6.0]
        sorted_group = [1.0,0.0,0.0,1.0,0.0,0.0]
        sorted_rank = array(sorted_rank)
        sorted_group = array(sorted_group)
        result = _compute_r_value(sorted_rank,sorted_group,4)
        self.assertEqual(result, .625)

    def test_anosim_p_test(self):
        """P-value should be .5 for this test"""
        group_list = {}
        samples, distmtx = parse_distmat(self.distmtx_txt)
        grouping, comment = parse_mapping_file_to_dict(self.mapping_txt)
        for sample in grouping:
            group_list[sample] = grouping[sample]["Treatment"]
            
        nrs = NonRandomShuffler()
            
        result, p_val = anosim_p_test(samples, distmtx,group_list, ntrials=3, randomfun=nrs.permutation)
        self.assertEqual(p_val, 0.5)
        
        
class NonRandomShuffler():
    """
    Class for testing the p_test function
        Since the p_test functino relies on randomness, it is necessary to use
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