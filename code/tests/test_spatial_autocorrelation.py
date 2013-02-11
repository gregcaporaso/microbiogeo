#!/usr/bin/env python
# File created on 8 Aug 2011

__author__ = "Andrew Cochran"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Andrew Cochran","Dan Knights"]
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "Dan Knights"
__email__ = "danknights@gmail.com"
__status__ = "Development"

from microbiogeo.spatial_autocorrelation import morans, _w_exponential,\
     _w_inverse, _s1, _s2, _s3, _s4, _s5, morans_variance
from cogent.util.unit_test import TestCase, main
from qiime.parse import parse_distmat, parse_mapping_file_to_dict
from numpy import array, zeros, max


class spatial_autocorrelationTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        self.distmtx_txt = """\tsam1\tsam2\tsam3\tsam4
        sam1\t0\t1\t5\t4
        sam2\t1\t0\t3\t2
        sam3\t5\t3\t0\t3
        sam4\t4\t2\t3\t0""".split('\n')
        
        self.w_inverse = """\tsam1\tsam2\tsam3\tsam4
        sam1\t0\t5\t1\t1.25
        sam2\t5\t0\t1.66666667\t2.5
        sam3\t1\t1.66666667\t0\t1.66666667
        sam4\t1.25\t2.5\t1.66666667\t0""".split('\n')

        self.w_exponential = """\tsam1\tsam2\tsam3\tsam4
        sam1\t0\t0.81873075\t0.36787944\t0.44932896
        sam2\t0.81873075\t0\t0.54881164\t0.67032005
        sam3\t0.36787944\t0.54881164\t0\t0.54881164
        sam4\t0.44932896\t0.67032005\t0.54881164\t0""".split('\n')
        self.mapping_txt="""#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tTreatment\tDOB\tDescription
        #Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
        sam1\tAGCACGAGCCTA\tYATGCTGCCTCCCGTAGGAGT\tControl\t8\tControl_mouse_I.D._354
        sam2\tAACTCGTCGATG\tYATGCTGCCTCCCGTAGGAGT\tControl\t3\tControl_mouse_I.D._355
        sam3\tACAGACCACTCA\tYATGCTGCCTCCCGTAGGAGT\tFast\t6\tControl_mouse_I.D._356
        sam4\tACCAGCGACTAG\tYATGCTGCCTCCCGTAGGAGT\tFast\t4\tControl_mouse_I.D._481""".split('\n')
        
        # Parse files
        samples, distmtx = parse_distmat(self.distmtx_txt)
        mapping_info, comment = parse_mapping_file_to_dict(self.mapping_txt)

        # Extract group column from mapping info
        self.x = zeros(len(samples))
        for i, sample in enumerate(samples):
            self.x[i] = mapping_info[sample]['DOB']
            
        # Normalize distances
        distmtx = distmtx/max(distmtx)

        # Weight w exponentialy
        self.w = _w_exponential(distmtx,1)
        self.sum_w = sum(sum(self.w))
        
    def test_morans(self):
        """Moran's should return -0.43078"""
        morans_res = morans(self.w,self.x)
        self.assertEqual(round(morans_res,5), -0.43078)
    
#    def test_varience(self):
#        """not yet implemented"""
#        i_value, v_value  = morans_variance(self.w, self.x)
#        self.assertEqual(round(v_value, 4) ,round(-18.8856654, 4))
        
    def test_w_inverse(self):
        """_w_inverse should return a matrix matching self.w_inverse"""
        samples, distmtx = parse_distmat(self.distmtx_txt)
        samples, w_inv = parse_distmat(self.w_inverse)
        distmtx = distmtx/max(distmtx)
        w = _w_inverse(distmtx)
        for i in range(len(w)):
            for j in range(len(w)):
                w[i][j] = round(w[i][j], 5)
                w_inv[i][j] = round(w_inv[i][j], 5)
        self.assertEqual(w,w_inv)

    def test_w_exponential(self):
        """_w_inverse should return a matrix matching self.w_exponential"""
        samples, distmtx = parse_distmat(self.distmtx_txt)
        samples, w_exp = parse_distmat(self.w_exponential)
        distmtx = distmtx/max(distmtx)
        w = _w_exponential(distmtx,1)
        for i in range(len(w)):
            for j in range(len(w)):
                w[i][j] = round(w[i][j], 5)
                w_exp[i][j] = round(w_exp[i][j], 5)
        self.assertEqual(w,w_exp)

    def test_s1(self):
        """Expected value: 8.23707"""
        s1 = _s1(self.w)
        self.assertEqual(round(s1, 4), round(8.23707, 4))

#    def test_s2(self):
#        """Expected value: 355.724909. Come back to this and do all the math?"""
#        s2 = _s2(self.w)
#        self.assertEqual(round(s2,4), round(355.724909,4))

    def test_s3(self):
        """Expected value: 1.573398"""
        s3 = _s3(self.w, [8,3,6,4], 5.25, 4)
        self.assertEqual(round(s3, 4), round(1.573398,4))

#    def test_s4(self):
#        """Expected value: -1226.203155026735"""
#        s4 = _s4(4, _s1(self.w), _s2(self.w), self.sum_w)
#        self.assertEqual(round(s4, 3), round(-1226.203155026735, 3))

    def test_s5(self):
        """Expected value: 220.414419465304"""
        s5 = _s5(_s1(self.w), 4, self.sum_w)
        self.assertEqual(round(s5, 4), round(220.414419465304, 4))
        
#run tests if called from command line
if __name__ == '__main__':
    main()
