#!/usr/bin/env python
# File created on 17 Mar 2011
# Modified by Logan Knecht on 03/13/2012 (Converted from script to Mantel class)
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso, Logan Knecht"]
__license__ = "GPL"
__version__ = "1.4.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"
 
from stats import CorrelationStats
from qiime.util import make_option
from qiime.parse import parse_distmat
from qiime.format import format_p_value_for_num_iters
from qiime.util import (parse_command_line_parameters, 
                        get_options_lookup,
                        make_compatible_distance_matrices)

from numpy import array, asarray, ravel, sqrt
from numpy.random import permutation

class Mantel(CorrelationStats):
    def __init__(self, samp_id_map, in_dm_fps, num_iterates):
        self.sample_id_map = samp_id_map
        self.input_dm_fps = in_dm_fps
        self.num_iterations = num_iterates
    
    def runAnalysis(self):
        resultsList = []
   
        for i,fp1 in enumerate(self.input_dm_fps):
            for fp2 in self.input_dm_fps[i+1:]:
                (dm1_labels, dm1), (dm2_labels, dm2) =\
                 make_compatible_distance_matrices(parse_distmat(open(fp1,'U')), parse_distmat(open(fp2,'U')), lookup=self.sample_id_map)
                if len(dm1_labels) < 2:
                    output_f.write('%s\t%s\t%d\tToo few samples\n' % (fp1,fp2,len(dm1_labels)))
                    continue

                m1, m2 = asarray(dm1), asarray(dm2)
                m1_flat = ravel(m1)
                size = len(m1)
                orig_stat = abs(self.pearson(m1_flat, ravel(m2)))
                better = 0
                for i in range(self.num_iterations):
                    p2 = self.permute_2d(m2, permutation(size))
                    r = abs(self.pearson(m1_flat, ravel(p2)))
                    if r >= orig_stat:
                        better += 1
                
                p = better

                p_str = format_p_value_for_num_iters(p,self.num_iterations)

                resultsList.append('%s\t%s\t%d\t%s\n' % (fp1,fp2,len(dm1_labels),p_str))
        return resultsList

    #This is a method was retrieved from the QIIME 1.4.0 release version, using amazon web services
    #Grabbed from the dir: /software/pycogent-1.5.1-release/lib/python2.7/site-packages/cogent/maths/stats
    #More specifically it was grabbed from the file called "test" 
    def pearson(self, x_items, y_items):
        """Returns Pearson correlation coefficient between x and y."""
        x_items, y_items = array(x_items), array(y_items)
        sum_x = sum(x_items)
        sum_y = sum(y_items)
        sum_x_sq = sum(x_items*x_items)
        sum_y_sq = sum(y_items*y_items)
        sum_xy = sum(x_items*y_items)
        n = len(x_items)
        try:
            r = 1.0 * ((n * sum_xy) - (sum_x * sum_y)) / \
               (sqrt((n * sum_x_sq)-(sum_x*sum_x))*sqrt((n*sum_y_sq)-(sum_y*sum_y)))
        except (ZeroDivisionError, ValueError, FloatingPointError): #no variation
            r = 0.0
        #check we didn't get a naughty value for r due to rounding error
        if r > 1.0:
            r = 1.0
        elif r < -1.0:
            r = -1.0
        return r
      
    #This is a method was retrieved from the QIIME 1.4.0 release version, using amazon web services
    #Grabbed from the dir: /software/pycogent-1.5.1-release/lib/python2.7/site-packages/cogent/maths/stats
    #More specifically it was grabbed from the file called "test"
    def permute_2d(self, m, p):
        """Performs 2D permutation of matrix m according to p."""
        return m[p][:, p]
