#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Test suite for the parse.py module."""

from cogent.util.unit_test import TestCase, main

from microbiogeo.parse import (parse_adonis_results,
                               parse_anosim_permanova_results,
                               parse_dbrda_results,
                               parse_mantel_results,
                               parse_morans_i_results,
                               parse_mrpp_results,
                               parse_partial_mantel_results,
                               parse_permdisp_results)

class ParseTests(TestCase):
    """Tests for the parse.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.anosim_results_str1 = anosim_results_str1.split('\n')
        self.anosim_results_str2 = anosim_results_str2.split('\n')

        self.adonis_results_str1 = adonis_results_str1.split('\n')
        self.adonis_results_str2 = adonis_results_str2.split('\n')

        self.mrpp_results_str1 = mrpp_results_str1.split('\n')

        self.dbrda_results_str1 = dbrda_results_str1.split('\n')

    def test_parse_anosim_permanova_results(self):
        """Test parsing anosim/permanova results file."""
        obs = parse_anosim_permanova_results(self.anosim_results_str1)
        self.assertFloatEqual(obs, (0.463253142506, 0.01))

        obs = parse_anosim_permanova_results(self.anosim_results_str2)
        self.assertFloatEqual(obs, (0.9375, None))

    def test_parse_adonis_results(self):
        """Test parsing adonis results file."""
        # Significant result.
        obs = parse_adonis_results(self.adonis_results_str1)
        self.assertFloatEqual(obs, (0.20273, 0.01))

        # Insignificant result.
        obs = parse_adonis_results(self.adonis_results_str2)
        self.assertFloatEqual(obs, (0.24408, 0.5))

    def test_parse_mrpp_results(self):
        """Test parsing mrpp results file."""
        obs = parse_mrpp_results(self.mrpp_results_str1)
        self.assertFloatEqual(obs, (0.07567, 0.01))

    def test_parse_dbrda_results(self):
        """Test parsing dbrda results file."""
        obs = parse_dbrda_results(self.dbrda_results_str1)
        print obs
        #self.assertFloatEqual(obs, (0.07567, 0.01))


anosim_results_str1 = """Method Name\tR-value\tP-value
ANOSIM\t0.463253142506\t0.01"""

anosim_results_str2 = """Method Name\tR-value\tP-value
ANOSIM\t0.9375\tToo few iters to compute p-value (num_iters=1)"""

adonis_results_str1 = """
Call:
adonis(formula = as.dist(qiime.data$distmat) ~ qiime.data$map[[opts$category]],      permutations = opts$num_permutations) 

                                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
qiime.data$map[[opts$category]]  6     6.635 1.10591  3.2632 0.20273   0.01 **
Residuals                       77    26.096 0.33891         0.79727          
Total                           83    32.731                 1.00000          
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1"""

adonis_results_str2 = """
Call:
adonis(formula = as.dist(qiime.data$distmat) ~ qiime.data$map[[opts$category]],      permutations = opts$num_permutations) 

                                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
qiime.data$map[[opts$category]]  1   0.44467 0.44467  2.2602 0.24408    0.5
Residuals                        7   1.37717 0.19674         0.75592       
Total                            8   1.82183                 1.00000       """

mrpp_results_str1 = """
Call:
mrpp(dat = as.dist(qiime.data$distmat), grouping = qiime.data$map[[opts$category]],      permutations = opts$num_permutations) 

Dissimilarity index: 
Weights for groups:  n 

Class means and counts:

      ENVO:forest ENVO:grassland ENVO:shrubland
delta 0.8133      0.8758         0.8456        
n     13          6              20            
      ENVO:Temperate broadleaf and mixed forest biome ENVO:Temperate grasslands
delta 0.7626                                          0.8554                   
n     19                                              11                       
      ENVO:Tropical and subtropical grasslands, savannas, and shrubland biome
delta 0.8392                                                                 
n     3                                                                      
      ENVO:Tropical humid forests
delta 0.7835                     
n     12                         

Chance corrected within-group agreement A: 0.07567 
Based on observed delta 0.8162 and expected delta 0.883 

Significance of delta: 0.01 
Based on  99  permutations"""

dbrda_results_str1 = """Analysis of Variance Table

Response: Distances
           Df Sum Sq   Mean Sq F value Pr(>F)
Groups      1 0.0119 0.0118769  2.0989 0.1479
Residuals 583 3.2990 0.0056587               

Permutation test for homogeneity of multivariate dispersions

No. of permutations: 999  

**** STRATA ****
Permutations are unstratified

**** SAMPLES ****
Permutation type: free 
Mirrored permutations for Samples?: No 

Response: Distances
           Df Sum Sq   Mean Sq      F N.Perm Pr(>F)
Groups      1 0.0119 0.0118769 2.0989    999  0.131
Residuals 583 3.2990 0.0056587                     

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
        female  male
female         0.132
male   0.14794  
"""


if __name__ == "__main__":
    main()
