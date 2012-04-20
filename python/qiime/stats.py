#!/usr/bin/env python
from __future__ import division

__author__ = "Michael Dwan"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Michael Dwan", "Logan Knecht", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Michael Dwan"
__email__ = "mdwan.tgen@gmail.com"
__status__ = "Development"

"""
This module provides functionality for the application of various statistical 
methods to QIIME formatted data sets.

The module provides classes, methods and functions that enable the user to
easily apply any number of statistical analyses and just as easily retrieve 
the results.
"""

from math import ceil, log, sqrt
from types import ListType

from cogent.cluster.metric_scaling import principal_coordinates_analysis
from cogent.maths.stats.test import pearson, permute_2d
from cogent.util.misc import combinate
from matplotlib import use
use('Agg', warn=False)
from matplotlib.pyplot import figure
from numpy import (add, array, asarray, asmatrix, dot, empty, finfo, matrix,
    newaxis, ravel, shape, square, std, transpose, zeros)
from numpy import min as np_min, max as np_max, sqrt as np_sqrt, sum as np_sum
from numpy.linalg import matrix_rank, qr, solve, svd
from numpy.random import permutation

from python.qiime.parse import DistanceMatrix, MetadataMap

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
        self.setData(mdmap, dm)
        self.setCategories(cats)
    
    def setData(self, new_map, new_dm):
        """Sets the instance's metadata map and distance matrix.

        Separate setter methods for the map and distance matrix are not
        provided because we need to be able to validate that the sample IDs
        match up between the two data structures.
      
        Arguments:
          new_map - A MetadataMap object instance.
          new_dm - A DistanceMatrix object instance.
        """
        if not isinstance(new_map, MetadataMap):
            raise TypeError('Invalid type: %s; not MetadataMap' %
                            new_map.__class__.__name__)
        if not isinstance(new_dm, DistanceMatrix):
            raise TypeError('Invalid type: %s; not DistanceMatrix' %
                            new_dm.__class__.__name__)
        if sorted(new_map.getSampleIds()) != sorted(new_dm.SampleIds):
            raise ValueError("The metadata map and distance matrix must have "
                "the same sample IDs.")
        self._metadata_map = new_map
        self.setDistanceMatrices([new_dm])

    def getMetadataMap(self):
        """Returns the instance's metadata map.
    
        The metadata map is returned as a MetadataMap class instance.
        """
        return self._metadata_map

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
            if not isinstance(el, str):
                raise TypeError("Invalid category: not of type 'string'")
            elif el not in self._metadata_map.getCategoryNames():
                raise ValueError("The category %s is not in the mapping file."
                    % el)
        self._categories = new_categories

    def getCategories(self):
        """Gets the instance's categories.
        
        Returns a list of mapping file category name strings.
        """
        return self._categories


class BioEnv(CategoryStats):
    """Class for BioEnv analysis."""

    def __init__(self, dm, metadata_map, cats):
        """Default constructor."""

        super(BioEnv, self).__init__(metadata_map, dm, cats)

    def runAnalysis(self):
        """Runs the BioEnv analysis on a distance matrix using specified
        metadata map categories.

        TODO: ADD COMMENTS

        Returns a dictionary which contains the resulting data. Keys:
            method_name - name of the statistical method
        
        Note: This code is loosely based on the implementation of BioEnv
        in the vegan package of R.
        """

        cats = self.getCategories()
        dm = self.getDistanceMatrices()[0]
        dm_flat = dm.flatten()

        row_count = dm.getSize()
        col_count = len(cats)
        sum = 0 
        stats = []
        for i in range(col_count+1):
            combo = list(combinate([j for j in range(0,col_count)], i))[1:]

            for c in range(len(combo)):
                cat_mat = self._make_cat_mat(cats, combo[c])
                cat_dm = self._derive_euclidean_dm(cat_mat, row_count)
                # stats.append(pearson(dm_flat, cat_dm.flatten()))
                stats.append(self._spearman_correlation(
                                 dm_flat, cat_dm.flatten()))

        sset = sorted(list(set(stats)))
        for s in sset:
            print s
        # print (2**col_count - 1)/2
       

    def _derive_euclidean_dm(self, cat_mat, dim):
        """Returns an n x n, euclidean distance matrix, where n = len(cats) """

        dm_labels = self.getDistanceMatrix().getSampleIds()
        res_mat = []
        for i in range(dim):
            res_mat.append([0 for k in range(dim)])
            for j in range(i):
                res_mat[i][j] = self._vector_dist(cat_mat[i], cat_mat[j])
                res_mat[j][i] = res_mat[i][j]

        return DistanceMatrix(asarray(res_mat), dm_labels, dm_labels)
    
    def _vector_dist(self, vec1, vec2):
        """Calculates the Euclidean distance between two vectors"""
        return sqrt(sum([(float(v1) - float(v2))**2 for v1,v2 in 
                            zip(vec1,vec2)]))


    def _make_cat_mat(self, cats, combo):
        """Returns a matrix with len(sample_ids) rows of columns pulled
        from category values, the number of columns for each category is
        determined by the current combination(combo)."""

        dm = self.getDistanceMatrix()
        md_map = self.getMetadataMap()
        res = []
        for i in combo:
            res.append(md_map.getCategoryValues(dm.getSampleIds(), cats[i]))

        return zip(*res)

    def _get_rank(self, data):
        """Ranks the elements of a list. Used in Spearman
        correlation
        """
        indices = range(len(data))
        ranks = range(1,len(data)+1)
        ranks.sort(key=lambda index:data[index-1])
        indices.sort(key=lambda index:data[index])
        data_len = len(data)
        i = 0
        while i < data_len: 
            j = i + 1
            val = data[indices[i]]
            while j < data_len and data[indices[j]] == val: 
                j += 1

            dup_ranks = j - i
            val = float(ranks[indices[i]]) + (dup_ranks-1)/2.0
            for k in range(i, i+dup_ranks):
                ranks[indices[k]] = val
            i += dup_ranks

        return ranks

    def _spearman_correlation(self, vec1, vec2):
        """Calculates the the Spearman distance of two vectors"""
        rank1 = self._get_rank(vec1)
        rank2 = self._get_rank(vec2)

        res = 0.0
        denom1 = 0.0
        denom2 = 0.0

        n = len(vec1)
        avg_rank = 0.5*(n-1)
        for i in range(n):
            res += rank1[i] * rank2[i]
            denom1 += rank1[i]**2
            denom2 += rank2[i]**2

        res = (res/n) - avg_rank**2
        denom1 = (denom1/n) - avg_rank**2
        denom2 = (denom2/n) - avg_rank**2

        if denom1 <= 0 or denom2 <= 0: 
            return 1.0

        res = 1.0 - (res/sqrt(denom1*denom2))

        return res


# if __name__ == '__main__':
#    dm = DistanceMatrix.parseDistanceMatrix(open('unweighted_unifrac_dm.txt'))
#    md_map = MetadataMap.parseMetadataMap(open('vars.txt'))
#    # dm = DistanceMatrix.parseDistanceMatrix(open('dm.txt'))
#    # md_map = MetadataMap.parseMetadataMap(open('vars2.txt'))


#    cats = ('TOT_ORG_CARB', 'SILT_CLAY', 'ELEVATION', 'SOIL_MOISTURE_DEFICIT', 'CARB_NITRO_RATIO', 'ANNUAL_SEASON_TEMP', 'ANNUAL_SEASON_PRECPT', 'PH', 'CMIN_RATE', 'LONGITUDE', 'LATITUDE')

#    bioenv = BioEnv(dm, md_map, cats)
#    # bioenv.runAnalysis()
#    a = (1,  2, 4, 3, 1, 6, 7, 8, 10, 4)
#    b = (2, 10, 20, 1, 3, 7, 5, 11, 6, 13)
#    print bioenv._spearman_correlation(a,b)




class DistanceBasedRda(CategoryStats):
    """Class for distance-based redundancy analysis."""

    def __init__(self, dm, metadata_map, category):
        """Default constructor."""
        if not isinstance(category, str):
            raise TypeError("The supplied category must be a string.")
        super(DistanceBasedRda, self).__init__(metadata_map, dm, [category])

    def getCategory(self):
        """Returns the single category of interest to this analysis."""
        return self.getCategories()[0]

    def setCategory(self, cat):
        """Sets the category of interest to this analysis.

        Arguments:
          cat - the category name (string). Must be present in the mapping
              file.
        """
        if not isinstance(cat, str):
            raise TypeError("The supplied category must be a string.")
        self.setCategories([cat])

    def runAnalysis(self):
        """Runs a distance-based redundancy analysis over the current data.

        TODO: add more useful comments

        Returns a dict containing the results. The following keys are set:
            method_name - name of the statistical method
        
        Note: This code is heavily based on the implementation of capscale/rda
        in R's vegan package.
        """
        results = {}
        results['method_name'] = "Distance-Based Redundancy Analysis"

        mdmap = self.getMetadataMap()
        dm = self.getDistanceMatrix()
        k = dm.getSize() - 1
        inertia_str = "Distance squared"

        if dm.getMax() >= 4:
            inertia_str = "mean " + inertia_str
            adjust = 1
        else:
            adjust = sqrt(k)

        points, eigs = principal_coordinates_analysis(
                asarray(dm.getDataMatrix()))

        # Order the axes in descending order based on eigenvalue. This code is
        # taken from QIIME's principal_coordinates.py library.
        idxs_descending = eigs.argsort()[::-1]
        points = points[idxs_descending]
        eigs = eigs[idxs_descending]

        # Drop the last dimension (axis) because it is invalid. PCoA is only
        # guaranteed to give up to (num_rows - 1) valid dimensions back, and we
        # get num_rows dimensions back. Since we've ordered the axes, we can
        # safely drop the last dimension.
        points = points[:-1]

        for eig in eigs:
            if eig < 0:
                raise ValueError("Encountered negative eigenvalue after "
                    "performing PCoA on the distance matrix. This might have "
                    "occurred if the distances are semi-metric or non-metric.")
        points = adjust * points
        if adjust == 1:
            eigs <- eigs / k

        group_membership = [mdmap.getCategoryValue(sid, self.getCategory()) \
                            for sid in dm.SampleIds]
        self._compute_rda(points, group_membership)

        return results

    def _compute_rda(self, point_matrix, group_membership):
        """Runs a traditional RDA analysis over the given data matrix.

        Returns all necessary information related to constrained/unconstrained
        axes (see vegan's rda method for more details).

        This method is heavily based on the RDA implementation found in R's
        vegan package.

        Arguments:
            point_matrix - numpy matrix where each row represents an axis and
                the columns represent points within that axis. Should be the
                output of principal_coordinates_analysis(). The columns will be
                in the order of the sample IDs that were in the original
                distance matrix.
            group_membership - a list of metadata mapping category values that
                indicate group membership for each of the samples. These values
                can be categorical or numerical and should be ordered the same
                as the sample IDs in the original distance matrix that was used
                as input to principal_coordinates_analysis(). Any strings
                representing categorical data will be converted to a number
                starting from zero (in the order they are encountered). Thus,
                ordinal categorical data must already be encoded as a number
                before being passed to this method. Otherwise, it will be
                interpretted as nominal categorical data.
        """
        zero = 0.0001
        cca = {}
        pcca = {}
        ca = {}

        # Transpose the point matrix so that it matches the type of matrix used
        # by R (i.e. rows are samples, and each column is a dimension/axis).
        point_matrix = asmatrix(point_matrix).T
        num_rows = point_matrix.shape[0] - 1
        points_bar = self._center_matrix(point_matrix)

        # Find the standard deviation of each column. Use df of 1 to match R's
        # sd function.
        points_bar_stdv = std(points_bar, axis=0, ddof=1)

        total_chi = sum(svd(points_bar, full_matrices=False,
                            compute_uv=False) ** 2) / num_rows

        # Do we need this?
        z_r = None

        factor = self._create_factor(group_membership)
        factor_r = self._center_matrix(factor)

        # Compute QR decomposition of the factor matrix.
        q, r = qr(factor_r)
        rank = matrix_rank(factor_r, tol=zero)
        qrank = rank
#        print "Q: "
#        print q
#        print "R: "
#        print r

        y = dot(q.T, points_bar)
        xQR = solve(r, y)
        #print xQR
        #print "points_bar: "
        #print points_bar
        #print "Q: "
        #print q
   #     Y <- qr.fitted(Q, Xbar)
   #     sol <- svd(Y)
   #     ## it can happen that rank < qrank
   #     rank <- min(rank, sum(sol$d > ZERO))
   #     sol$d <- sol$d/sqrt(NR)
   #     ax.names <- paste("RDA", 1:length(sol$d), sep = "")
   #     colnames(sol$u) <- ax.names
   #     colnames(sol$v) <- ax.names
   #     names(sol$d) <- ax.names
   #     rownames(sol$u) <- rownames(X)
   #     rownames(sol$v) <- colnames(X)


    def _center_matrix(self, mat):
        """Returns a column-centered version of the matrix.

        A column-centered matrix will have the mean of each column in the
        input matrix subtracted from each element in that column.

        Returns a numpy matrix of type float.

        This code is partially derived from http://stackoverflow.com/a/8917508.

        Arguments:
            mat - should be a numpy matrix
        """
        # Convert to float type because testing this with an integer numpy
        # matrix gave a resulting matrix of ints, which is probably not what we
        # want here.
        centered = mat.copy().astype(float)
        centered -= centered.sum(0) / centered.shape[0]
        return centered

    def _create_factor(self, group_membership):
        """Transforms group membership list into a factor.
        
        The factor is basically a column vector. Any categorical data (strings)
        are transformed into a numeric representation. This does not respect
        ordinal categorical data; it will simply represent the first category
        it encounters as a zero. This code is mimicing R's 'factor' data
        structure.

        The return value is a numpy column vector (i.e. nx1 matrix), where n is
        len(group_membership). For example, say we have the following input:
            ["Fast", "Fast", "Control", "Fast", "Control"]

        The output would be (written as a python list, but it will be a column
        vector in actuality):
            [0, 0, 1, 0, 1]

        The returned column vector will contain floats.

        Arguments:
            group_membership - a list of categorical data values (strings) or a
                numeric list if it represents numeric data. If the list is
                numeric, it will simply be returned as a numpy column vector
                (i.e. no transformation).
        """
        try:
            factor = map(float, group_membership)
        except:
            # We have categorical data represented as strings.
            factor = []
            category_symbol = 0
            conversion = {}
            for ele in group_membership:
                if ele not in conversion:
                    conversion[ele] = category_symbol
                    category_symbol += 1
                factor.append(conversion[ele])
        return matrix(factor).T


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

        This method overrides its parent. Only two distance matrices may be
        set, and the smallest allowable size is 3x3 (for Pearson correlation).
        
        Arguments:
            matrices - list of exactly two DistanceMatrix objects.
        """
        if len(matrices) != 2:
            raise ValueError("Can only set exactly two distance matrices for "
                             "a Mantel correlogram analysis.")
        if matrices[0].getSize() < 3 or matrices[1].getSize() < 3:
            raise ValueError("Both distance matrices must be at least 3x3.")
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

        # Find the number of lower triangular elements (excluding the
        # diagonal).
        num_dists = dm_size * (dm_size - 1) // 2

        # Use Sturge's rule to determine the number of distance classes.
        num_classes = int(ceil(1 + log(num_dists, 2)))

        # Create the matrix of distance classes. Each element in the matrix
        # contains what distance class the original element is in. Also find
        # the distance class indices, which are the midpoints in each distance
        # class.
        dist_class_matrix, class_indices = self._find_distance_classes(geo_dm,
            num_classes)
        
        # Start assembling the results.
        results = {}
        results['method_name'] = 'Mantel Correlogram'
        results['class_index'] = []
        results['num_dist'] = []
        results['mantel_r'] = []
        results['mantel_p'] = []

        # Create a model matrix for each distance class, then compute a Mantel
        # test using it and the original eco distance matrix. A model matrix
        # contains ones for each element that is in the current distance class,
        # and zeros otherwise (zeros on the diagonal as well).
        for class_num in range(num_classes):
            results['class_index'].append(class_indices[class_num])
            model_matrix = zeros([dm_size, dm_size], dtype=int)
            for i in range(dm_size):
                for j in range(dm_size):
                    curr_ele = dist_class_matrix[i][j]
                    if curr_ele == class_num and i != j:
                        model_matrix[i][j] = 1
            model_matrix = DistanceMatrix(model_matrix, geo_dm.SampleIds,
                                          geo_dm.SampleIds)

            # Count the number of distances in the current distance class.
            num_distances = int(model_matrix.sum())
            results['num_dist'].append(num_distances)
            if num_distances == 0:
                results['mantel_r'].append(None)
                results['mantel_p'].append(None)
            else:
                row_sums = model_matrix.sum(axis='observation')
                row_sums = map(int, row_sums)
                has_zero_sum = 0 in row_sums

                # Only stop running Mantel tests if we've gone through half of
                # the distance classes and at least one row has a sum of zero
                # (i.e. the sample doesn't have any distances that fall in the
                # current class).
                if not (class_num > ((num_classes // 2) - 1) and has_zero_sum):
                    mantel_test = Mantel(model_matrix, eco_dm,
                            self.getNumPermutations(), tail_type='greater')
                    mantel_test_results = mantel_test.runAnalysis()
                    p_val, orig_stat, perm_stats = (
                            mantel_test_results['p_value'],
                            mantel_test_results['r_value'],
                            mantel_test_results['perm_stats'])
                    results['mantel_r'].append(-orig_stat)

                    # The mantel function produces a one-tailed p-value
                    # (H1: r>0). Here, compute a one-tailed p-value in the
                    # direction of the sign.
                    if orig_stat < 0:
                        perm_sum = sum([1 for ps in perm_stats \
                            if ps <= orig_stat]) + 1
                        p_val = perm_sum / (self.getNumPermutations() + 1)
                    results['mantel_p'].append(p_val)
                else:
                    results['mantel_r'].append(None)
                    results['mantel_p'].append(None)

        # Correct p-values for multiple testing using Bonferroni correction
        # (non-progressive).
        num_tests = len([p_val for p_val in results['mantel_p'] \
                         if p_val is not None])
        corrected_p_vals = [min(p * num_tests, 1) \
                            for p in results['mantel_p'][0:num_tests]]
        corrected_p_vals.extend([None] * (num_classes - num_tests))
        results['mantel_p_corr'] = corrected_p_vals

        # Construct a correlogram of distance class versus mantel correlation
        # statistic and fill in each point that is statistically significant.
        results['correlogram_plot'] = self._generate_correlogram(
            results['class_index'], results['mantel_r'],
            results['mantel_p_corr'])
        return results

    def _find_distance_classes(self, dm, num_classes):
        """Computes a distance class matrix and distance class midpoints.
        
        Returns a matrix of the same dimensions as the input matrix but each
        element indicates which distance class (0..num_classes-1) the original
        element belongs to. The diagonal will always have a value of -1,
        indicating that it is not apart of any distance class. Also returns a
        list of distance class midpoints.

        Distance classes are determined by the minimum and maximum values in
        the input matrix and the number of specified classes.

        Arguments:
            dm - the input DistanceMatrix object to compute distance classes on
            num_classes - the number of desired distance classes
        """

        if num_classes < 1:
            raise ValueError("Cannot have fewer than one distance class.")

        # Compute the breakpoints of the distance classes based on the number
        # of specified classes and the ranges of values in the lower triangular
        # portion of the distance matrix (excluding the diagonal).
        dm_lower_flat = dm.flatten()
        break_points = self._find_break_points(np_min(dm_lower_flat),
            np_max(dm_lower_flat), num_classes)

        # Find the class indices (the midpoints between breakpoints).
        class_indices = []
        for bp_index, break_point in enumerate(break_points[0:num_classes]):
            next_bp = break_points[bp_index + 1]
            class_indices.append(break_point + (0.5 * (next_bp - break_point)))

        # Create the matrix of distance classes. Every element in the matrix
        # tells what distance class the original element belongs to.
        size = dm.getSize()
        dist_class_matrix = empty([size, size], dtype=int)
        for i in range(size):
            for j in range(size):
                if i != j:
                    curr_ele = dm[i][j]
                    bps = [(k - 1) for k, bp in enumerate(break_points) \
                        if bp >= curr_ele]
                    dist_class_matrix[i][j] = min(bps)
                else:
                    dist_class_matrix[i][j] = -1
        return dist_class_matrix, class_indices

    def _find_break_points(self, start, end, num_classes):
        """Finds the points to break a range into equal width classes.

        Returns a list of floats indicating breakpoints in the range.

        Arguments:
            start - the minimum value in the range
            end - the maximum value in the range
            num_classes - the number of classes to break the range into
        """
        if start >= end:
            raise ValueError("Cannot find breakpoints because the starting "
                "point is greater than or equal to the ending point.")
        if num_classes < 1:
            raise ValueError("Cannot have fewer than one distance class.")

        width = (end - start) / num_classes
        break_points = [start + width * class_num \
            for class_num in range(num_classes)]
        break_points.append(float(end))

        # Move the first breakpoint a little bit to the left. Machine epsilon
        # is take from:
        # http://en.wikipedia.org/wiki/Machine_epsilon#
        #     Approximation_using_Python
        epsilon = finfo(float).eps
        break_points[0] = break_points[0] - epsilon
        return break_points

    def _generate_correlogram(self, class_indices, mantel_stats,
            corrected_p_vals):
        """Generates a matplotlib plot of the Mantel correlogram.

        Returns a matplotlib Figure instance, which can then be manipulated
        further or saved to a file as necessary.

        Arguments:
            class_indices - list of distance class indices (for the x-axis)
            mantel_stats - list of Mantel r stats (for the y-axis)
            corrected_p_vals - list of corrected p-values (for filling in
                points to indicate significance)
        """
        # Plot distance class index versus mantel correlation statistic.
        fig = figure()
        ax = fig.add_subplot(111)
        ax.plot(class_indices, mantel_stats, 'ks-', mfc='white', mew=1)

        # Fill in each point that is significant (based on alpha).
        signif_classes = []
        signif_stats = []
        for idx, p_val in enumerate(corrected_p_vals):
            if p_val <= self.getAlpha():
                signif_classes.append(class_indices[idx])
                signif_stats.append(mantel_stats[idx])
        ax.plot(signif_classes, signif_stats, 'ks', mfc='k')

        ax.set_title("Mantel Correlogram")
        ax.set_xlabel("Distance class index")
        ax.set_ylabel("Mantel correlation statistic")
        return fig


class Mantel(CorrelationStats):
    """Class for the Mantel matrix correlation statistical method.

    This class provides the functionality to run a Mantel analysis on two
    distance matrices.
    
    TODO: Put plain english explanation here, have Damien explain it.
    """
    def __init__(self, dm1, dm2, permutations, tail_type="two sided"):
        """Constructs a new Mantel instance.

        Arguments:
            dm1 - first DistanceMatrix object to be compared
            dm2 - second DistanceMatrix object to be compared
            permutations - the number of times to permute the distance matrix
                while calculating the p-value
            tail_type - the type of Mantel test to perform (i.e. hypothesis
                test). Can be "two sided", "less", or "greater"
        """
        super(Mantel, self).__init__([dm1, dm2])
        self.setNumPermutations(permutations)
        self.setTailType(tail_type)

    def runAnalysis(self):
        results = self._mantel_test()

        resultsDict = {}
        resultsDict['method_name'] = "Mantel"
        resultsDict['dm1'] = self.getDistanceMatrices()[0]
        resultsDict['dm2'] = self.getDistanceMatrices()[1]
        resultsDict['num_perms'] = self.getNumPermutations() 
        resultsDict['p_value'] = results[0]
        resultsDict['r_value'] = results[1]
        resultsDict['perm_stats'] = results[2]
        resultsDict['tail_type'] = self.getTailType()

        return resultsDict

    def _mantel_test(self):
        """Runs a Mantel test on the current distance matrices.
    
        Returns the p-value, Mantel correlation statistic, and a list of Mantel
        correlation statistics for each permutation test. The currently set
        tail type and number of permutations will be used to run the test.

        Note: this method was taken from the development version of PyCogent as
        we needed access to different tail types and the currently released
        version of PyCogent does not support this. Once this functionality is
        available in the version of PyCogent supported by QIIME, we should
        remove this method and use the one in PyCogent instead. This method
        isn't exactly the same as the PyCogent implementation because it has
        been adapted to use the class members and DistanceMatrix objects, but
        in essence it is the same implementation.
        """
        m1, m2 = self.getDistanceMatrices()
        n = self.getNumPermutations()
        alt = self.getTailType()

        # Get a flattened list of lower-triangular matrix elements (excluding
        # the diagonal) in column-major order. Use these values to calculate
        # the correlation statistic.
        m1_flat, m2_flat = m1.flatten(True), m2.flatten(True)
        orig_stat = pearson(m1_flat, m2_flat)

        # Run our permutation tests so we can calculate a p-value for the test.
        better = 0
        perm_stats = []
        for i in range(n):
            m1_perm_data = permute_2d(m1, permutation(m1.getSize()))
            m1_perm = DistanceMatrix(m1_perm_data, m1.getSampleIds(),
                m1.getSampleIds())
            m1_perm_flat = m1_perm.flatten()
            r = pearson(m1_perm_flat, m2_flat)

            if alt == 'two sided':
                if abs(r) >= abs(orig_stat):
                    better += 1
            else:
                if ((alt == 'greater' and r >= orig_stat) or
                    (alt == 'less' and r <= orig_stat)):
                    better += 1
            perm_stats.append(r)
        return (better + 1) / (n + 1), orig_stat, perm_stats

    def getTailType(self):
        """Returns the tail type being used for the Mantel test."""
        return self._tail_type

    def setTailType(self, tail_type):
        """Sets the tail type that will be used for the Mantel test.
        
        Valid types are 'two sided', 'less', or 'greater'.
        """
        if tail_type not in ("two sided", "greater", "less"):
            raise ValueError("Unrecognized alternative hypothesis (tail "
                             "type). Must be either 'two sided', 'greater', "
                             "or 'less'.")
        self._tail_type = tail_type

    def getNumPermutations(self):
        """Returns the number of permutations used."""
        return self._num_perms

    def setNumPermutations(self, num_perms):
        """Sets the number of permutations to be used."""
        if num_perms >= 0:
            self._num_perms = num_perms
        else:
            raise ValueError("The number of permutations cannot be a negative "
                             "value.")

    def setDistanceMatrices(self, matrices):
        """Sets the distance matrices being used for the analysis.
        
        Only two distance matrices are allowed in a Mantel test.
        """
        if len(matrices) != 2:
            raise ValueError("Can only set exactly two distance matrices for "
                             "a Mantel test.")
        super(Mantel, self).setDistanceMatrices(matrices)


class PartialMantel(CorrelationStats):
    def __init__(self, dm1, dm2, cdm, num_perms):
        super(PartialMantel, self).__init__([dm1, dm2, cdm])
        self.setNumPermutations(num_perms)

    def getNumPermutations(self):
        return self._num_perms

    def setNumPermutations(self, num_perms):
        if num_perms >= 0:
            self._num_perms = num_perms
        else:
            raise ValueError("The number of permutations cannot be a negative value.")

    def setDistanceMatrices(self, matrices):
        if len(matrices) != 3:
            raise ValueError("partial Mantel analysis requires two distance matrices and a third control matrix.")
        super(PartialMantel, self).setDistanceMatrices(matrices)

    def runAnalysis(self):
        """Run a partial Mantel test on the current distance matrices and control matrix.
        
        Credit: The code herein is based loosely on the implementation found in the Vegan
                package of the R language and software libraries.
        """

        # calculate the correlation statistic. 
        corr = lambda rxy, rxz, ryz: (rxy - rxz*ryz)/(sqrt(1 - rxz**2)*sqrt(1 - ryz**2))

        # Load initial/placeholder values in the results dictionary.
        res = {}
        res['method_name'] = 'partial Mantel'
        res['mantel_r'] = None
        res['mantel_p'] = None

        perm_num = self.getNumPermutations()

        dm1 = self.getDistanceMatrices()[0]
        dm2 = self.getDistanceMatrices()[1]
        cdm = self.getDistanceMatrices()[2]
        dm_sizes = dm1.getSize()

        dm1_flat = dm1.flatten()
        dm2_flat = dm2.flatten()
        cdm_flat = cdm.flatten()
        
        # Get the initial r-values before permuting.
        rval1 = pearson(dm1_flat, dm2_flat)
        rval2 = pearson(dm1_flat, cdm_flat)
        rval3 = pearson(dm2_flat, cdm_flat)

        # Calculate the orginal test statistic (p-value)
        orig_stat = corr(rval1, rval2, rval3)

        # Calculate permuted r-values and p-values, storing
        # them for use in the calculation of the final statistic.
        perm_stats = [0 for i in range(perm_num)]
        numerator = 0
        for i in range(0,perm_num):
            # Permute the first distance matrix and calculate new
            # r and p-values 
            p1 = permute_2d(dm1, permutation(dm1.getSize()))
            dm1_perm = DistanceMatrix(p1, dm1.SampleIds, dm1.SampleIds)
            dm1_perm_flat = dm1_perm.flatten()
            rval1 = pearson(dm1_perm_flat, dm2_flat)
            rval2 = pearson(dm1_perm_flat, cdm_flat)
            perm_stats.append(corr(rval1, rval2, rval3))

            # Sum the permuted statistics for calculation of the final statistic.
            if perm_stats[-1] >= orig_stat:
              numerator += perm_stats[-1]

        
        # Load the final statistics into the result dictionary.
        res['mantel_r'] = orig_stat
        res['mantel_p'] = (numerator + 1) / (perm_num + 1)
        
        return res
