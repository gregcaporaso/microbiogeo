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
        self.adonis_results_str1 = adonis_results_str1.split('\n')
        self.adonis_results_str2 = adonis_results_str2.split('\n')

    def test_parse_adonis_results(self):
        """Test parsing adonis results file."""
        # Significant result.
        obs = parse_adonis_results(self.adonis_results_str1)
        self.assertFloatEqual(obs, (0.20273, 0.01))

        # Insignificant result.
        obs = parse_adonis_results(self.adonis_results_str2)
        self.assertFloatEqual(obs, (0.24408, 0.5))


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


if __name__ == "__main__":
    main()
