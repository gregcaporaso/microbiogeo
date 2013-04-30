#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Test suite for the method.py module."""

from cogent.util.unit_test import TestCase, main

from microbiogeo.method import (AbstractStatMethod, Adonis, Anosim, Best,
                                Dbrda, Mantel, MantelCorrelogram, MoransI,
                                Mrpp, PartialMantel, Permanova, Permdisp,
                                QiimeStatMethod, UnparsableFileError,
                                UnparsableLineError)

class AbstractStatMethodTests(TestCase):
    """Tests for the AbstractStatMethod class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.inst = AbstractStatMethod()

    def test_parse(self):
        """Test raises error."""
        self.assertRaises(NotImplementedError, self.inst.parse, 'foo')

    def test_parse_float(self):
        """Test parsing float strings."""
        obs = self.inst.parse_float('0.045')
        self.assertFloatEqual(obs, 0.045)

        obs = self.inst.parse_float('1.0')
        self.assertFloatEqual(obs, 1.0)

        obs = self.inst.parse_float('0')
        self.assertFloatEqual(obs, 0.0)

        obs = self.inst.parse_float('42')
        self.assertFloatEqual(obs, 42.0)

        obs = self.inst.parse_float('5', 0)
        self.assertFloatEqual(obs, 5.0)

        obs = self.inst.parse_float('-9', max_val=10)
        self.assertFloatEqual(obs, -9.0)

    def test_parse_float_invalid_input(self):
        """Test parsing invalid float strings raises errors."""
        self.assertRaises(ValueError, self.inst.parse_float, 'foo')
        self.assertRaises(ValueError, self.inst.parse_float, '5.0', 0, 1)
        self.assertRaises(ValueError, self.inst.parse_float, '-0.1', 0, 1)
        self.assertRaises(ValueError, self.inst.parse_float, '-0.1', 0)
        self.assertRaises(ValueError, self.inst.parse_float, '11.2',
                          max_val=10)


class QiimeStatMethodTests(TestCase):
    """Tests for the QiimeStatMethod class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.inst = QiimeStatMethod()

        self.anosim_results_str1 = anosim_results_str1.split('\n')
        self.anosim_results_str2 = anosim_results_str2.split('\n')

    def test_parse(self):
        """Test parsing QIIME stats results file."""
        obs = self.inst.parse(self.anosim_results_str1)
        self.assertFloatEqual(obs, (0.463253142506, 0.01))

        self.assertRaises(UnparsableLineError, self.inst.parse,
                          self.anosim_results_str2)


class AnosimTests(TestCase):
    """Nothing to test."""
    pass


class PermanovaTests(TestCase):
    """Nothing to test."""
    pass


class AdonisTests(TestCase):
    """Tests for the Adonis class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.inst = Adonis()

        self.adonis_results_str1 = adonis_results_str1.split('\n')
        self.adonis_results_str2 = adonis_results_str2.split('\n')
        self.adonis_results_str3 = adonis_results_str3.split('\n')

    def test_parse(self):
        """Test parsing adonis results file."""
        # Significant result.
        obs = self.inst.parse(self.adonis_results_str1)
        self.assertFloatEqual(obs, (0.20273, 0.01))

        # Insignificant result.
        obs = self.inst.parse(self.adonis_results_str2)
        self.assertFloatEqual(obs, (0.24408, 0.5))

        # Invalid format.
        self.assertRaises(UnparsableFileError, self.inst.parse,
                          self.adonis_results_str3)


class MrppTests(TestCase):
    """Tests for the Mrpp class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.inst = Mrpp()

        self.mrpp_results_str1 = mrpp_results_str1.split('\n')
        self.mrpp_results_str2 = mrpp_results_str2.split('\n')

    def test_parse(self):
        """Test parsing mrpp results file."""
        obs = self.inst.parse(self.mrpp_results_str1)
        self.assertFloatEqual(obs, (0.07567, 0.01))

        # Invalid format.
        self.assertRaises(UnparsableFileError, self.inst.parse,
                          self.mrpp_results_str2)


class DbrdaTests(TestCase):
    """Tests for the Dbrda class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.inst = Dbrda()

        self.dbrda_results_str1 = dbrda_results_str1.split('\n')

    def test_parse(self):
        """Test parsing dbrda results file."""
        obs = self.inst.parse(self.dbrda_results_str1)
        self.assertFloatEqual(obs, (0.2786, 0.010101))


class PermdispTests(TestCase):
    """Tests for the Permdisp class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.inst = Permdisp()

        self.permdisp_results_str1 = permdisp_results_str1.split('\n')

    def test_parse(self):
        """Test parsing permdisp results file."""
        obs = self.inst.parse(self.permdisp_results_str1)
        self.assertFloatEqual(obs, (2.0989, 0.131))


class MantelTests(TestCase):
    """Tests for the Mantel class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.inst = Mantel()

        self.mantel_results_str1 = mantel_results_str1.split('\n')

    def test_parse(self):
        """Test parsing mantel results file."""
        obs = self.inst.parse(self.mantel_results_str1)
        self.assertFloatEqual(obs, (1.0, 0.01))


class PartialMantelTests(TestCase):
    """Tests for the PartialMantel class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.inst = PartialMantel()

        self.partial_mantel_results_str1 = \
                partial_mantel_results_str1.split('\n')

    def test_parse(self):
        """Test parsing partial mantel results file."""
        obs = self.inst.parse(self.partial_mantel_results_str1)
        self.assertFloatEqual(obs, (0.5, 0.01))


class MantelCorrelogramTests(TestCase):
    """Tests for the MantelCorrelogram class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.inst = MantelCorrelogram()

    def test_parse(self):
        """Test raises error."""
        self.assertRaises(NotImplementedError, self.inst.parse, 'foo')


class MoransITests(TestCase):
    """Tests for the MoransI class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.inst = MoransI()

        self.morans_i_results_str1 = morans_i_results_str1.split('\n')

    def test_parse(self):
        """Test parsing moran's i results file."""
        obs = self.inst.parse(self.morans_i_results_str1)
        self.assertFloatEqual(obs, (-0.06005486, 4.442088e-05))


class BestTests(TestCase):
    """Tests for the Best class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.inst = Best()

    def test_parse(self):
        """Test raises error."""
        self.assertRaises(NotImplementedError, self.inst.parse, 'foo')


anosim_results_str1 = """Method name\tR statistic\tp-value\tNumber of permutations
ANOSIM\t0.463253142506\t0.01\t99"""

anosim_results_str2 = """Method name\tR statistic\tp-value\tNumber of permutations
ANOSIM\t0.9375\tToo few iters to compute p-value (num_iters=1)\t1"""

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

adonis_results_str3 = """
Call:
adonis(formula = as.dist(qiime.data$distmat) ~ qiime.data$map[[opts$category]],      permutations = opts$num_permutations) 

                                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
foo.data$map[[opts$category]]  1   0.44467 0.44467  2.2602 0.24408    0.5
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

mrpp_results_str2 = """
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

Significance of foo: 0.01 
Based on  99  permutations"""

dbrda_results_str1 = """Call: capscale(formula = as.dist(qiime.data$distmat) ~ factor, data =
factors.frame)

               Inertia Proportion Rank
Total         159.1762                
Real Total    165.4413     1.0000     
Constrained    46.0873     0.2786   19
Unconstrained 119.3540     0.7214  371
Imaginary      -6.2651             213
Inertia is squared Unknown distance 

Eigenvalues for constrained axes:
    CAP1     CAP2     CAP3     CAP4     CAP5     CAP6     CAP7     CAP8 
14.72239 10.95891  8.89776  3.26489  2.89957  1.41151  0.87627  0.69475 
    CAP9    CAP10    CAP11    CAP12    CAP13    CAP14    CAP15    CAP16 
 0.40960  0.35446  0.29999  0.24395  0.20137  0.18342  0.17567  0.15110 
   CAP17    CAP18    CAP19 
 0.13347  0.11498  0.09327 

Eigenvalues for unconstrained axes:
  MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
12.480  5.688  4.495  3.722  3.331  2.814  2.279  2.153 
(Showed only 8 of all 371 unconstrained eigenvalues)


Permutation test for capscale 

Call: capscale(formula = as.dist(qiime.data$distmat) ~ factor, data =
factors.frame)
Permutation test for all constrained eigenvalues
Pseudo-F:\t 11.48258 (with 19, 565 Degrees of Freedom)
Significance:\t 0.010101 
Based on 98 permutations under reduced model.
"""

permdisp_results_str1 = """Analysis of Variance Table

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

mantel_results_str1 = """# Number of entries refers to the number of rows (or cols) retained in each
# distance matrix after filtering the distance matrices to include only those
# samples that were in both distance matrices. p-value contains the correct
# number of significant digits.
DM1\tDM2\tNumber of entries\tMantel r statistic\tp-value\tNumber of permutations\tTail type
/Users/jrideout/analysis/overview_tutorial/wf_bdiv_even146/unweighted_unifrac_dm.txt\t/Users/jrideout/analysis/overview_tutorial/wf_bdiv_even146/unweighted_unifrac_dm.txt\t9\t1.00000\t0.01\t100\ttwo sided"""

partial_mantel_results_str1 = """# Number of entries refers to the number of rows (or cols) retained in each
# distance matrix after filtering the distance matrices to include only those
# samples that were in both distance matrices. p-value contains the correct
# number of significant digits.
DM1\tDM2\tCDM\tNumber of entries\tMantel r statistic\tp-value\tNumber of permutations\tTail type
/Users/jrideout/analysis/overview_tutorial/wf_bdiv_even146/unweighted_unifrac_dm.txt\t/Users/jrideout/analysis/overview_tutorial/wf_bdiv_even146/unweighted_unifrac_dm.txt\t/Users/jrideout/analysis/overview_tutorial/wf_bdiv_even146/unweighted_unifrac_dm.txt\t9\t0.50000\t0.01\t100\tgreater"""

morans_i_results_str1 = """$observed
[1] -0.06005486

$expected
[1] -0.125

$sd
[1] 0.01590547

$p.value
[1] 4.442088e-05
"""


if __name__ == "__main__":
    main()
