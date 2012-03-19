#!/usr/bin/env python
from __future__ import division

__author__ = "Michael Dwan"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Michael Dwan", "Logan Knecht", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Michael Dwan"
__email__ = "mgd25@nau.edu"
__status__ = "Development"

"""
This module provides functionality for the application of various statistical 
methods to QIIME formatted data sets.

The module provides classes, methods and functions that enable the user to
easily apply any number of statistical analyses and easily retrieve the 
results.
"""

from cogent.maths.stats.test import pearson, permute_2d
from math import ceil, log
from matplotlib import use
use('Agg', warn=False)
from matplotlib.pyplot import figure
from numpy import array, asarray, empty, ravel
from numpy.random import permutation
from types import ListType

from python.qiime.parse import DistanceMatrix

class GradientStats(object):
    """Top-level, abstract base class for gradient statistical analyses.
    
    It exists to allow for extension to non-matrix based analyses.
    """

    def __init__(self):
        pass

    def runAnalysis(self):
        raise NotImplementedError("Method not implemented by abstract base.")


class DistanceMatrixStats(GradientStats):
    """Base class for distance matrix-based statistical methods.
    
    It is the parent class of CorrelationStats and CategoryStats. They extend
    from this class in order to provided consistent method use for similar
    classes. More specifically to provide the same runAnalysis method
    throughout.
    """

    def __init__(self, distmats):
        """Default constructor.

        Arguments:
          distmats - a list of DistanceMatrix objects
        """
        super(DistanceMatrixStats, self).__init__()
        self.setDistanceMatrices(distmats)

    def getDistanceMatrices(self):
        """Returns the list of distance matrices."""
        return self._distmats
  
    def setDistanceMatrices(self, matrices):
        """Sets the list of distance matrices to the supplied list.

        Arguments:
          matrices - the new list of distance matrices being assigned
        """
        if not isinstance(matrices, ListType):
            raise TypeError("The item passed in as the new list was not a "
                            "list data type.")
        self._distmats = matrices


class CorrelationStats(DistanceMatrixStats):
    """Base class for distance matrix correlation statistical methods.
    
    It is subclassed by correlation methods like Partial Mantel and Mantel that
    compare two or more distance matrices. A valid instance of CorrelationStats
    must have at least one distance matrix, and all distance matrices must have
    matching dimensions and sample IDs (i.e. matching row/column labels).
    """

    def __init__(self, distmats):
        """Default constructor.
        
        Arguments:
          distmats - a list of DistanceMatrix objects
        """
        super(CorrelationStats, self).__init__(distmats)

    def setDistanceMatrices(self, matrices):
        if len(matrices) < 1:
            raise ValueError("Must provide at least one distance matrix.")

        size = matrices[0].getSize()
        sample_ids = matrices[0].SampleIds
        for dm in matrices:
            if dm.getSize() != size:
                raise ValueError("All distance matrices must have the same "
                                 "number of rows and columns.")
            if dm.SampleIds != sample_ids:
                raise ValueError("All distance matrices must have matching "
                                 "sample IDs.")
        super(CorrelationStats, self).setDistanceMatrices(matrices)


class CategoryStats(DistanceMatrixStats):
    """Base class for categorical statistical analyses."""

    def __init__(self, mdmap, dm, cats):
        """Default constructor."""
        super(CategoryStats, self).__init__([dm])
        self._metadata_map = mdmap
        self._categories = cats
    
    def setMetadataMap(self, new_map):
        """Sets the instance's metadata map to a new MetadataMap instance.
      
        Arguments:
          new_map - A MetadataMap object instance.
        """
        if not isinstance(new_map, self.__class__):
            raise TypeError('Invalid type: %s; not MetadataMap' %
                            new_map.__class__.__name__)
        self._metadata_map = new_map

    def getMetadataMap(self):
        """Returns the instance's metadata map.
    
        The metadata map is returned as a MetadataMap class instance.
        """
        return self._metadata_map

    def setDistanceMatrix(self, new_distmat):
        """Sets the instance's distance matrix.
    
        Arguments:
          new_distmat - A DistanceMatrix object instance.
        """
        if not isinstance(new_distmat, self.__class__):
            raise TypeError('Invalid type: %s; not DistanceMatrix' %
                            new_distmat.__class__.__name__)
        self.setDistanceMatrices([new_distmat])

    def getDistanceMatrix(self):
        """Gets the instance's distance matrix.
    
        The distance matrix is returned as a DistanceMatrix class instance.
        """
        return self.getDistanceMatrices()[0]

    def setCategories(self, new_categories):
        """Sets the instance's list of categories to a new list of strings
        representing categories in a QIIME mapping file.
    
        Arguments:
          new_categories - A list of category name strings.
        """
        for el in new_categories:
          if not isinstance(el, self.__class__):
            raise TypeError('Invalid category: not of type "string"')
        self._categories = new_categories

    def getCategories(self):
        """Gets the instance's categories.
        
        Returns a list of mapping file category name strings.
        """
        return self._categories


class MantelCorrelogram(CorrelationStats):
    """Class for the Mantel correlogram statistical method.

    This class provides the functionality to run a Mantel correlogram analysis
    on two distance matrices. In a nutshell, the distances are split into
    distance classes and a Mantel test is run over each distance class. A
    Mantel correlogram is created, which is basically a plot of distance
    classes versus Mantel statistics.

    Uses Sturge's rule to determine the number of distance classes, and
    Pearson's method to compute the correlation at each distance class. The
    corrected p-values are computed using Bonferroni correction.
    """

    def __init__(self, eco_dm, geo_dm, num_perms, alpha=0.05):
        """Constructs a new MantelCorrelogram instance.

        Arguments:
            eco_dm - a DistanceMatrix object representing the ecological
                distances between samples (e.g. UniFrac distance matrix).
            geo_dm - a DistanceMatrix object representing some other distance
                measure between samples (most commonly geographical distances,
                but could also be distances in pH, temperature, etc.).
            num_perms - the number of permutations to use when computing the
                p-values.
            alpha - the alpha value to use when marking the Mantel
                correlogram plot for significance.
        """
        super(MantelCorrelogram, self).__init__([eco_dm, geo_dm])
        self.setNumPermutations(num_perms)
        self.setAlpha(alpha)

    def getNumPermutations(self):
        """Returns the number of permutations to use in the Mantel tests."""
        return self._num_perms

    def setNumPermutations(self, num_perms):
        """Sets the number of permutations to use in the Mantel tests.
        
        Arguments:
            num_perms - the number of permutations. This value must be greater
                than or equal to zero.
        """
        if num_perms >= 0:
            self._num_perms = num_perms
        else:
            raise ValueError("The number of permutations cannot be negative.")

    def getAlpha(self):
        """Returns the alpha value."""
        return self._alpha

    def setAlpha(self, alpha):
        """Sets the alpha value.

        Arguments:
            alpha - the value of alpha. Must be between 0 and 1, inclusive.
        """
        if alpha >= 0 and alpha <= 1:
            self._alpha = alpha
        else:
            raise ValueError("Alpha must be between 0 and 1.")

    def setDistanceMatrices(self, matrices):
        """Sets the distance matrices to use in the Mantel correlogram test.

        This method overrides its parent.
        
        Arguments:
            matrices - list of exactly two DistanceMatrix objects.
        """
        if len(matrices) != 2:
            raise ValueError("Can only set exactly two distance matrices for "
                             "a Mantel correlogram analysis.")
        super(MantelCorrelogram, self).setDistanceMatrices(matrices)

    def runAnalysis(self):
        """Runs a Mantel correlogram test over the current distance matrices.

        Returns a dict containing the results. The following keys are set:
            method_name - name of the statistical method
            class_index - list of distance class indices (the center of each
                distance class)
            num_dist - list of the number of distances in each distance class
            mantel_r - list of the Mantel r statistics for each distance class
            mantel_p - list of the p-values for each distance class
            mantel_p_corr - list of the p-values for each distance class,
                corrected for multiple tests
            correlogram_plot - a matplotlib Figure object containing the
                correlogram
        
        Note: This code is heavily based on the implementation of
        mantel.correlog in R's vegan package.
        """
        eco_dm = self.getDistanceMatrices()[0]
        geo_dm = self.getDistanceMatrices()[1]
        dm_size = eco_dm.getSize()

        # Find the number of lower/upper triangular elements (discounting the
        # diagonal).
        num_dists = dm_size * (dm_size - 1) // 2

        # Use Sturge's rule to determine the number of distance classes.
        num_classes = int(ceil(1 + log(num_dists, 2)))

        # Compute the breakpoints based on the number of distance classes.
        flattened_lower = geo_dm.flatten()
        start_point = min(flattened_lower)
        end_point = max(flattened_lower)
        width = (end_point - start_point) / num_classes
        break_points = []
        for class_num in range(num_classes):
            break_points.append(start_point + width * class_num)
        break_points.append(end_point)

        # Move the first breakpoint a little bit to the left. Machine epsilon
        # is take from:
        # http://en.wikipedia.org/wiki/Machine_epsilon#
        #     Approximation_using_Python
        epsilon = 2.2204460492503131e-16
        break_points[0] = break_points[0] - epsilon

        # Find the class indices (the midpoints between breakpoints).
        class_indices = []
        for bp_index, break_point in enumerate(break_points[0:num_classes]):
            next_bp = break_points[bp_index + 1]
            class_index = break_point + (0.5 * (next_bp - break_point))
            class_indices.append(class_index)

        # Create the matrix of distance classes.
        flattened_geo_dm = geo_dm.flatten(lower=False)
        dist_class_matrix = []
        for dm_ele in flattened_geo_dm:
            bps = [i for i, bp in enumerate(break_points) if bp >= dm_ele]
            dist_class_matrix.append(min(bps))

        # Start assembling the vectors of results.
        class_index = [None]
        num_dist = [None]
        mantel_r = [None]
        mantel_p = [None]

        # Create a model-matrix for each distance class, then compute a Mantel
        # test.
        for class_idx in range(num_classes):
            class_index.append(class_indices[class_idx])
            model_matrix_list = [0] * (dm_size ** 2)
            for ele_idx, ele in enumerate(dist_class_matrix):
                # Fix this hack so we don't need to add 1.
                if ele == (class_idx + 1):
                    model_matrix_list[ele_idx] = 1
            model_matrix = empty([dm_size, dm_size], dtype=int)
            # Convert vector into matrix, setting diagonal to zero.
            for idx, ele in enumerate(model_matrix_list):
                col_num = idx // dm_size
                row_num = idx % dm_size
                if row_num == col_num:
                    model_matrix[row_num][col_num] = 0
                else:
                    model_matrix[row_num][col_num] = ele
            model_matrix = DistanceMatrix(model_matrix, geo_dm.SampleIds,
                                          geo_dm.SampleIds)
            num_distances = int(model_matrix.sum())
            num_dist.append(num_distances)
            if num_distances == 0:
                mantel_r.append(None)
                mantel_p.append(None)
            else:
                row_sums = model_matrix.sum(axis='observation')
                row_sums = map(int, row_sums)
                # Fix this hack so we don't need to add 1.
                has_zero_sum = False
                for row_sum in row_sums:
                    if row_sum == 0:
                        has_zero_sum = True
                        break
                # Only stop running Mantel tests if we've gone through half of
                # the distance classes and at least one row has a sum of zero
                # (i.e. the sample doesn't have any distances that fall in the
                # current class).
                if not ((class_idx + 1) > (num_classes // 2) and has_zero_sum):
                    temp_p_val, orig_stat, perm_stats = self._mantel(
                        model_matrix._data, eco_dm._data,
                        self.getNumPermutations())
                    mantel_r.append(-orig_stat)

                    # The mantel() function produces a one-tailed p-value
                    # (H1: r>0). Here, compute a one-tailed p-value in the
                    # direction of the sign.
                    if orig_stat < 0:
                        perm_sum = sum([1 for ps in perm_stats \
                            if ps <= orig_stat]) + 1
                        temp_p_val = perm_sum / (self.getNumPermutations() + 1)
                    mantel_p.append(temp_p_val)
                else:
                    mantel_r.append(None)
                    mantel_p.append(None)

        results = {}
        results['method_name'] = 'Mantel Correlogram'
        results['class_index'] = class_index[1:]
        results['num_dist'] = num_dist[1:]
        results['mantel_r'] = mantel_r[1:]
        results['mantel_p'] = mantel_p[1:]
        
        # List mantel_p starts with a None value.
        mantel_p = mantel_p[1:]
        num_tests = len([p_val for p_val in mantel_p if p_val is not None])

        # Correct p-values for multiple testing using Bonferroni correction
        # (non-progressive).
        corrected_p_vals = [min(p * num_tests, 1) \
                            for p in mantel_p[0:num_tests]]
        corrected_p_vals.extend([None] * (num_classes - num_tests))
        results['mantel_p_corr'] = corrected_p_vals

        # Construct a matplotlib plot of distance class versus mantel
        # correlation statistic.
        fig = figure()
        ax = fig.add_subplot(111)
        ax.plot(results['class_index'], results['mantel_r'], 'ks-',
                mfc='white', mew=1)
        # Fill in each point that is significant (corrected p-value <= 0.05).
        signif_classes = []
        signif_stats = []
        for idx, p_val in enumerate(results['mantel_p_corr']):
            if p_val <= self.getAlpha():
                signif_classes.append(results['class_index'][idx])
                signif_stats.append(results['mantel_r'][idx])
        ax.plot(signif_classes, signif_stats, 'ks', mfc='k')

        ax.set_title("Mantel Correlogram")
        ax.set_xlabel("Distance class index")
        ax.set_ylabel("Mantel correlation statistic")
        results['correlogram_plot'] = fig

        #fig.savefig('mantel_correlogram.png', format='png')

        return results
            
    def _mantel(self, dm1, dm2, num_perms):
        """Runs a Mantel test over the supplied distance matrices.

        Returns a tuple containing the p-value, Mantel r statistic, and the
        Mantel r statistic for each permutation.

        The first distance matrix is the one that is permuted when calculating
        the p-value. The p-value is based on a one-tailed test (H1: r>0). The
        Mantel r statistic is computed using Pearson's correlation method.

        This code is based on R's vegan::mantel function.

        Arguments:
            dm1 - DistanceMatrix object
            dm2 - DistanceMatrix object
            num_perms - the number of permutations, must be >= 0
        """
        dm1, dm2 = asarray(dm1), asarray(dm2)
        samp_ids = self.getDistanceMatrices()[0].SampleIds
        dm1_dm = DistanceMatrix(dm1, samp_ids, samp_ids)
        dm2_dm = DistanceMatrix(dm2, samp_ids, samp_ids)
        dm1_flat = dm1_dm.flatten()
        dm2_flat = dm2_dm.flatten()
        size = dm1_dm.getSize()
        orig_stat = pearson(dm1_flat, dm2_flat)
        better = 0
        perm_stats = []
        for i in range(num_perms):
            p1 = permute_2d(dm1, permutation(size))
            p1_dm = DistanceMatrix(p1, samp_ids, samp_ids)
            p1_flat = p1_dm.flatten()
            r = pearson(p1_flat, dm2_flat)
            perm_stats.append(r)
            if r >= orig_stat:
                better += 1
        return (better + 1) / (num_perms + 1), orig_stat, perm_stats

#    def _mantel2(self, m1, m2, n):
#        """Compares two distance matrices. Reports P-value for correlation."""
#        m1, m2 = asarray(m1), asarray(m2)
#        m1_flat = ravel(m1)
#        size = len(m1)
#        #orig_stat = abs(pearson(m1_flat, ravel(m2)))
#        orig_stat = pearson(m1_flat, ravel(m2))
#        better = 0
#        perm_stats = []
#        for i in range(n):
#            #p2 = m2[permutation(size)][:, permutation(size)]
#            p2 = permute_2d(m2, permutation(size))
#            #r = abs(pearson(m1_flat, ravel(p2)))
#            r = pearson(m1_flat, ravel(p2))
#            perm_stats.append(r)
#            if r >= orig_stat:
#                better += 1
#        return better/n, orig_stat, perm_stats
