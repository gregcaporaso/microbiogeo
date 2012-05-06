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

"""
This module provides functionality for the application of various statistical
methods to QIIME-formatted datasets.

The module provides an API that allows users to easily apply any number of
statistical analyses and just as easily retrieve the results. The module also
provides a hierarchy of statistical classes that can be inherited from to
create new statistical method implementations.
"""

from types import ListType

from cogent.cluster.metric_scaling import principal_coordinates_analysis
from cogent.maths.stats.test import pearson, permute_2d
from cogent.util.misc import combinate
from matplotlib import use
use('Agg', warn=False)
from matplotlib.pyplot import figure
from numpy import (arange, argsort, array, asarray, asmatrix, ceil, dot, empty,
                   finfo, log2, matrix, mean, ones, random, sqrt, std, tri,
                   unique, zeros)
from numpy import min as np_min, max as np_max
from numpy.linalg import matrix_rank, qr, solve, svd
from numpy.random import permutation

from qiime.stats import CategoryStats
from qiime.util import DistanceMatrix, MetadataMap


class BioEnv(CategoryStats):
    """Class for BioEnv analysis."""

    def __init__(self, dm, metadata_map, cats):
        """Default constructor."""

        super(BioEnv, self).__init__(metadata_map, [dm], cats, num_dms=1)

    def __call__(self):
        """Runs the BioEnv analysis on a distance matrix using specified
        metadata map categories.

        TODO: ADD COMMENTS

        Returns a dictionary which contains the resulting data. Keys:
            method_name - name of the statistical method
        """

        res = super(BioEnv, self).__call__()
        cats = self.Categories
        dm = self.DistanceMatrices[0]
        dm_flat = dm.flatten()
        dm_flat_ranked = self._get_rank(dm_flat)

        row_count = dm.Size
        col_count = len(cats)
        sum = 0
        stats = [(-777777777, '') for c in range(col_count+1)]
        for i in range(1, col_count+1):
            combo = list(combinate([j for j in range(1,col_count+1)], i))

            for c in range(len(combo)):
                cat_mat = self._make_cat_mat(cats, combo[c])
                cat_dm = self._derive_euclidean_dm(cat_mat, row_count)
                cat_dm_flat_ranked = self._get_rank(cat_dm.flatten())
                r = self._spearman_correlation(dm_flat_ranked,
                                               cat_dm_flat_ranked, ranked=True)
                if r > stats[i-1][0]:
                    stats[i-1] = (r, ','.join(str(s) for s in combo[c]))

        # for s in stats[1:]:
        #     print s
        # print (2**col_count - 1)/2

        res['method_name'] = 'BioEnv'
        res['num_vars'] = col_count
        res['vars'] = ['%s = %d' % (name,val+1) for val,name in enumerate(cats)]
        res['bioenv_rho_vals'] = stats[:-1]

        return res


    def _derive_euclidean_dm(self, cat_mat, dim):
        """Returns an n x n, euclidean distance matrix, where n = len(cats) """

        dm_labels = self.DistanceMatrices[0].SampleIds
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

        dm = self.DistanceMatrices[0]
        md_map = self.MetadataMap
        res = []
        for i in combo:
            res.append(md_map.getCategoryValues(dm.SampleIds, cats[i-1]))

        return zip(*res)

    def _get_rank(self, data):
        """Ranks the elements of a list. Used in Spearman
        correlation
        """
        indices = range(len(data))
        ranks = range(1,len(data)+1)
        indices.sort(key=lambda index:data[index])
        ranks.sort(key=lambda index:indices[index-1])
        data_len = len(data)
        i = 0
        ties = 0
        while i < data_len:
            j = i + 1
            val = data[indices[i]]
            try:
                val += 0
            except TypeError:
                raise(TypeError)

            while j < data_len and data[indices[j]] == val:
                j += 1

            dup_ranks = j - i
            val = float(ranks[indices[i]]) + (dup_ranks-1)/2.0
            for k in range(i, i+dup_ranks):
                ranks[indices[k]] = val
            i += dup_ranks

            ties += dup_ranks-1

        return ranks, ties

    def _spearman_correlation(self, vec1, vec2, ranked=False):
        """Calculates the the Spearman distance of two vectors"""
        try:
            temp = len(vec1)
        except ValueError:
            raise(ValueError, 'First input vector is not a list.')

        try:
            temp = len(vec2)
        except ValueError:
            raise(ValueError, 'Second input vector is not a list.')

        if len(vec1) == 0 or len(vec2) == 0:
            raise(ValueError, 'One or both input vectors has/have zero elements')

        if len(vec1) != len(vec2):
            raise(ValueError, 'Vector lengths must be equal')

        if not ranked:
            rank1, ties1 = self._get_rank(vec1)
            rank2, ties2 = self._get_rank(vec2)
        else:
            rank1, ties1 = vec1
            rank2, ties2 = vec2

        if ties1 == 0 and ties2 == 0:
            n = len(rank1)
            sum_sqr = sum([(x-y)**2 for x,y in zip(rank1,rank2)])
            rho = 1 - (6*sum_sqr/(n*(n**2 - 1)))
            return rho

        avg = lambda x: sum(x)/len(x)

        x_bar = avg(rank1)
        y_bar = avg(rank2)

        numerator = sum([(x-x_bar)*(y-y_bar) for x,y in zip(rank1, rank2)])
        denominator = sqrt(sum([(x-x_bar)**2 for x in rank1])*sum([(y-y_bar)**2 for y in rank2]))

        rho = numerator/denominator
        return rho

if __name__ == '__main__':
    dm = DistanceMatrix.parseDistanceMatrix(open('unweighted_unifrac_dm.txt'))
    md_map = MetadataMap.parseMetadataMap(open('map.txt'))
   # dm = DistanceMatrix.parseDistanceMatrix(open('dm.txt'))
   # md_map = MetadataMap.parseMetadataMap(open('vars2.txt'))


    cats = ['TOT_ORG_CARB', 'SILT_CLAY', 'ELEVATION', 'SOIL_MOISTURE_DEFICIT', 'CARB_NITRO_RATIO', 'ANNUAL_SEASON_TEMP', 'ANNUAL_SEASON_PRECPT', 'PH', 'CMIN_RATE', 'LONGITUDE', 'LATITUDE']
    # cats = ['LONGITUDE', 'LATITUDE']

    bioenv = BioEnv(dm, md_map, cats)
    print bioenv()
#     a = [1,2,4,3,1,6,7,8,10,4]
#     b = [2,10,20,1,3,7,5,11,6,13]
#     x = (1,  2, 4, 3, 1, 6, 7, 8, 10, 4, 100, 2, 3, 77)
#     y = (2, 10, 20, 1, 3, 7, 5, 11, 6, 13, 5, 6, 99, 101)
#     r = (1,2,4,5,2,2,4,3,1,4)
#     s = (2,3,5,4,2,2,3,4,3,2)
#     u = (1,2,3,4,5,6,7,8,9)
#     v = (10, 11, 4, 2, 9, 33, 1, 5, 88)
    # print bioenv._spearman_correlation(a,b)
    # print bioenv._spearman_correlation(x,y)
    # print bioenv._spearman_correlation(r,s)
    # print bioenv._spearman_correlation(u,v)

   # x = c(1,  2, 4, 3, 1, 6, 7, 8, 10, 4, 100, 2, 3, 77)
   # y = c(2, 10, 20, 1, 3, 7, 5, 11, 6, 13, 5, 6, 99, 101)



class DistanceBasedRda(CategoryStats):
    """Class for distance-based redundancy analysis."""

    def __init__(self, dm, metadata_map, category):
        """Default constructor."""
        if not isinstance(category, str):
            raise TypeError("The supplied category must be a string.")
        super(DistanceBasedRda, self).__init__(metadata_map, [dm], [category],
                                               num_dms=1)

    def __call__(self, num_perms=999):
        """Runs a distance-based redundancy analysis over the current data.

        TODO: add more useful comments

        Returns a dict containing the results. The following keys are set:
            method_name - name of the statistical method

        Note: This code is heavily based on the implementation of capscale/rda
        in R's vegan package.
        """
        results = super(DistanceBasedRda, self).__call__()
        results['method_name'] = "Distance-Based Redundancy Analysis"

        mdmap = self.MetadataMap
        dm = self.DistanceMatrices[0]
        k = dm.Size - 1
        inertia_str = "Distance squared"

        if dm.max() >= 4:
            inertia_str = "mean " + inertia_str
            adjust = 1
        else:
            adjust = sqrt(k)

        points, eigs = principal_coordinates_analysis(dm.DataMatrix)

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

        group_membership = [mdmap.getCategoryValue(sid, self.Categories[0]) \
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
