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

from math import ceil, log, sqrt
from types import ListType

from cogent.cluster.metric_scaling import principal_coordinates_analysis
from cogent.maths.stats.test import pearson, permute_2d
from cogent.util.misc import combinate
from matplotlib import use
use('Agg', warn=False)
from matplotlib.pyplot import figure
from numpy import (arange, argsort, array, asarray, asmatrix, dot, empty,
                   finfo, matrix, mean, std, tri, zeros)
from numpy import min as np_min, max as np_max
from numpy.linalg import matrix_rank, qr, solve, svd
from numpy.random import permutation

from python.qiime.parse import DistanceMatrix, MetadataMap


class DistanceMatrixStats(object):
    """Base class for distance matrix-based statistical methods.

    This class provides an interface to setting and accessing an arbitrary
    number of distance matrices. Users of this class can optionally specify the
    number of allowable distance matrices and their minimum allowable size (the
    default is no restrictions on either of these).

    It is the parent class of CorrelationStats and CategoryStats.
    """

    def __init__(self, distmats, num_dms=-1, min_dm_size=-1):
        """Default constructor.
        
        Initializes an instance with the provided list of distance matrices.

        Arguments:
            distmats - a list of DistanceMatrix objects
            num_dms - the exact number of allowable distance matrices. If -1
                (the default), there is no restriction on how many distance
                matrices the user can set
            min_dm_size - the minimum size that all distance matrices must have
                that are stored by this instance. If -1, no size restriction
        """
        self._num_dms = num_dms
        self._min_dm_size = min_dm_size
        self.setDistanceMatrices(distmats)

    def getDistanceMatrices(self):
        """Returns the list of distance matrices."""
        return self._distmats
  
    def setDistanceMatrices(self, distmats):
        """Sets the list of distance matrices to the supplied list.

        Arguments:
            distmats - the new list of distance matrices being assigned
        """
        if not isinstance(distmats, ListType):
            raise TypeError("The item passed in as the new list was not a "
                            "list data type.")
        if self._num_dms >= 0 and len(distmats) != self._num_dms:
            raise ValueError("Cannot set %d distance matrices. Must provide "
                             "exactly %d distance matrices." % (len(distmats),
                             self._num_dms))
        for dm in distmats:
            if not isinstance(dm, DistanceMatrix):
                raise TypeError('Invalid type: %s; expected DistanceMatrix' %
                                dm.__class__.__name__)
            elif self._min_dm_size >= 0 and dm.getSize() < self._min_dm_size:
                raise ValueError("Distance matrix of size %dx%d is smaller "
                                 "than the minimum allowable distance matrix "
                                 "size of %dx%d for this analysis." %
                                 (dm.getSize(), dm.getSize(),
                                  self._min_dm_size, self._min_dm_size))
        self._distmats = distmats

    def runAnalysis(self):
        """Runs the statistical method and returns relevant results.

        The return value of this method is a python dictionary with arbitrary
        key/value pairs of results, since each statistical method returns
        different results.

        This method is not implemented and should be implemented by subclasses.
        It is a good idea to call the parent class' runAnalysis() method first
        to obtain any results from the parent and then add more results to the
        dict that is obtained from the parent.
        """
        raise NotImplementedError("Method not implemented by abstract base.")


class PermutationStats(object):
    """Base class for permutation-based statistical methods.

    This generic class should be inherited from if the statistical method is
    permutation-based (e.g. Mantel permutes a number of times to compute a
    p-value).
    """

    def __init__(self, num_perms):
        """Default constructor creates an instance with the number of perms.

        Arguments:
            num_perms - the number of permutations. This value must be greater
                than or equal to zero
        """
        self.setNumPermutations(num_perms)

    def getNumPermutations(self):
        """Returns the number of permutations to use."""
        return self._num_perms

    def setNumPermutations(self, num_perms):
        """Sets the number of permutations to use.

        Arguments:
            num_perms - the number of permutations. This value must be greater
                than or equal to zero
        """
        try :
            num_perms = int(num_perms)
        except:
            raise TypeError("The number of permutations supplied must be an "
                            "integer.")
        if num_perms >= 0:
            self._num_perms = num_perms
        else:
            raise ValueError("The number of permutations cannot be negative.")

    def runAnalysis(self):
        """Runs the statistical method and returns relevant results.

        The return value of this method is a python dictionary with arbitrary
        key/value pairs of results, since each statistical method returns
        different results.

        This method is not implemented and should be implemented by subclasses.
        It is a good idea to call the parent class' runAnalysis() method first
        to obtain any results from the parent and then add more results to the
        dict that is obtained from the parent.
        """
        raise NotImplementedError("Method not implemented by abstract base.")


class CorrelationStats(DistanceMatrixStats):
    """Base class for distance matrix correlation statistical methods.

    It is subclassed by correlation methods such as partial Mantel and Mantel
    that compare two or more distance matrices.

    A valid instance of CorrelationStats must have at least one distance
    matrix, and all distance matrices must have matching dimensions and sample
    IDs (i.e. matching row/column labels). This check is in place to prevent
    the accidental comparison on two distance matrices that have sample IDs in
    different orders. Essentially, all of the distance matrices must be
    "compatible".

    Users of this class can optionally specify the number of allowable distance
    matrices and their minimum allowable size (the default is no restrictions
    on either of these).
    """

    def __init__(self, distmats, num_dms=-1, min_dm_size=-1):
        """Default constructor.

        Creates a new instance with the provided list of distance matrices.

        Arguments:
            distmats - a list of DistanceMatrix objects
            num_dms - the exact number of allowable distance matrices. If -1
                (the default), there is no restriction on how many distance
                matrices the user can set
            min_dm_size - the minimum size that all distance matrices must have
                that are stored by this instance. If -1, no size restriction
        """
        super(CorrelationStats, self).__init__(distmats, num_dms, min_dm_size)

    def setDistanceMatrices(self, distmats):
        """Sets the list of distance matrices to the supplied list.

        This method overrides the parent method and enforces more checks to
        ensure that at least one distance matrix is provided and that all of
        the distance matrices are compatible.

        Arguments:
            distmats - the new list of distance matrices being assigned
        """
        super(CorrelationStats, self).setDistanceMatrices(distmats)
        if len(distmats) < 1:
            raise ValueError("Must provide at least one distance matrix.")

        size = distmats[0].getSize()
        sample_ids = distmats[0].SampleIds
        for dm in distmats:
            if dm.getSize() != size:
                raise ValueError("All distance matrices must have the same "
                                 "number of rows and columns.")
            if dm.SampleIds != sample_ids:
                raise ValueError("All distance matrices must have matching "
                                 "sample IDs.")


class CategoryStats(object):
    """Base class for categorical statistical analyses.

    It is subclassed by categorical statistical methods such as DB-RDA or BEST.
    Categorical statistical methods usually have some categorical grouping of
    samples, and the significance of this grouping is usually what is tested.
    For example, are treatment samples significantly different from control
    samples? This is not always the case (e.g. DB-RDA is an ordination
    technique), but most of the categorical methods follow this general design.

    A valid instance of CategoryStats must have one distance matrix and an
    accompanying metadata map.
    """

    def __init__(self, mdmap, dm, cats):
        """Default constructor.

        Creates a new instance with the provided distance matrix, metadata map,
        and list of categories.

        Arguments:
            mdmap - a MetadataMap instance
            dm - a DistanceMatrix instance
            cats - a list of strings denoting categories in the metadata map
                that will be used by this analysis (i.e. the grouping
                variable(s))
        """
        self.setData(mdmap, dm)
        self.setCategories(cats)

    def setData(self, new_mdmap, new_dm):
        """Sets the instance's metadata map and distance matrix.

        Separate setter methods for the map and distance matrix are not
        provided because we need to be able to validate that the sample IDs
        match up between the two data structures. It seems to be a safe
        assumption that if the user is changing the distance matrix, the
        metadata map also should be changed at the same time.

        Arguments:
            new_mdmap - A MetadataMap object instance
            new_dm - A DistanceMatrix object instance
        """
        if not isinstance(new_mdmap, MetadataMap):
            raise TypeError('Invalid type: %s; not MetadataMap' %
                            new_mdmap.__class__.__name__)
        if not isinstance(new_dm, DistanceMatrix):
            raise TypeError('Invalid type: %s; not DistanceMatrix' %
                            new_dm.__class__.__name__)
        if not self.compatibleSampleIds(new_dm, new_mdmap):
            raise ValueError("The metadata map and distance matrix must have "
                             "the same sample IDs.")
        self._metadata_map = new_mdmap
        self._dm = new_dm

    def getMetadataMap(self):
        """Returns the instance's metadata map.
 
        The metadata map is returned as a MetadataMap class instance.
        """
        return self._metadata_map

    def getDistanceMatrix(self):
        """Gets the instance's distance matrix.

        The distance matrix is returned as a DistanceMatrix class instance.
        """
        return self._dm

    def setCategories(self, new_categories):
        """Sets the instance's list of categories.

        Arguments:
            new_categories - A list of category name strings. These must be
                present in the current metadata map
        """
        if not isinstance(new_categories, ListType):
            raise TypeError("The supplied categories must be a list of "
                            "strings.")
        for new_cat in new_categories:
            if not isinstance(new_cat, str):
                raise TypeError("Invalid category: not of type 'string'")
            elif new_cat not in self._metadata_map.getCategoryNames():
                raise ValueError("The category %s is not in the mapping file."
                    % new_cat)
        self._categories = new_categories

    def getCategories(self):
        """Gets the instance's categories.

        Returns a list of mapping file category name strings.
        """
        return self._categories

    def compatibleSampleIds(self, dm, mdmap):
        """Returns True if the sample IDs are the same in both structures.

        This method will return True if the sample IDs in the distance matrix
        are exactly the same as the sample IDs in the metadata map. Ordering of
        samples is not taken into account.

        Arguments:
            dm - the DistanceMatrix instance to compare
            mdmap - the MetadataMap instance to compare
        """
        same = False
        if sorted(dm.getSampleIds()) == sorted(mdmap.getSampleIds()):
            same = True
        return same

    def runAnalysis(self):
        """Runs the statistical method and returns relevant results.

        The return value of this method is a python dictionary with arbitrary
        key/value pairs of results, since each statistical method returns
        different results.

        This method is not implemented and should be implemented by subclasses.
        It is a good idea to call the parent class' runAnalysis() method first
        to obtain any results from the parent and then add more results to the
        dict that is obtained from the parent.
        """
        raise NotImplementedError("Method not implemented by abstract base.")


class Anosim(CategoryStats, PermutationStats):
    """Class for ANOSIM categorical statistical analysis.

    This code is heavily based on Andrew Cochran's original procedural version.
    """

    def __init__(self, mdmap, dm, cat, num_perms):
        CategoryStats.__init__(self, mdmap, dm, [cat])
        PermutationStats.__init__(self, num_perms)

    def runAnalysis(self):
        ntrials = self.getNumPermutations()
        category = self.getCategories()[0]
        samples = self.getDistanceMatrix().getSampleIds()
        distmtx = self.getDistanceMatrix().getDataMatrix()

        # Array to hold value of p-tests
        r_value_permunations = zeros(ntrials)

        group_hash = {}
        for samp_id in samples:
            group_hash[samp_id] = self.getMetadataMap().getCategoryValue(
                    samp_id, category)

        # Calculate the R value
        r_value = self._anosim(samples, distmtx, group_hash)

        # Main loop to run the p-tests
        for i in xrange(ntrials):
            # Randomize grouping
            grouping_random = []
            for sample in samples:
                grouping_random.append(group_hash[sample])
            grouping_random = permutation(grouping_random)

            for j, sample in enumerate(samples):
                group_hash[sample] = grouping_random[j]
            r_value_permunations[i] = self._anosim(samples, distmtx,
                                                   group_hash)

        # Calculate the p-value and return
        p_value = (sum(r_value_permunations >= r_value) + 1) / (ntrials + 1)

        return {'method_name': 'ANOSIM', 'r_value': r_value,
                'p_value': p_value}

    def _anosim(self, samples, distmtx, group_hash):
        """Computes ANOSIM on the current data, returning the R value.

        The R value is between -1 and 1 and indicates the strength of the
        grouping.
        """
        distmtx = asarray(distmtx)
        # Local Variable
        matrix_size = len(distmtx)

        # Create grouping matrix
        with_between = zeros(distmtx.shape)
        for i, i_sample in enumerate(samples):
            group_list_i_sample = group_hash[i_sample]
            for j, j_sample in enumerate(samples):
                if group_list_i_sample == group_hash[j_sample]:
                    with_between[i][j] = 1
        
        # Extract upper triangle
        differences = distmtx[tri(len(distmtx)) == 0]
        grouping = with_between[tri(len(with_between)) == 0]
        
        # Sort extracted data
        sorted_differences = []
        sorted_grouping = []
        for idx in argsort(differences):
             sorted_differences.append(differences[idx])
             sorted_grouping.append(grouping[idx])
        sorted_differences = array(sorted_differences)
        sorted_grouping = array(sorted_grouping)
        
        # Account for rank ties, then compute r value
        rank_list = arange(1,len(sorted_differences) + 1)
        adjusted_rank_list = self._remove_ties(sorted_differences, rank_list)
        result = self._compute_r_value(adjusted_rank_list, sorted_grouping,
                                       matrix_size)
        return result

    def _remove_ties(self, sorted_diffs, ranks):
        """
        Replaces repeat values with the average of them
        
            PARAMETERS
            sorted_diffs: array of the sorted differences
            ranks: array containing the ranks of each of the differences
            
            RETURNS
            result: array of the adjusted ranks
        """
        # Local Variables
        result = []
        tie_list = []
        tie_count = 0
        tie_flag = 0
        
        # Main Loop
        for i in range(len(sorted_diffs)-1):
            # Store state information
            curr = sorted_diffs[i]
            next = sorted_diffs[i+1]
            rank_val = ranks[i]
            
            # A tie has not occured yet
            if tie_flag == 0:
                # Check for a tie
                if curr == next:
                    tie_count = tie_count + 1
                    tie_list.append(rank_val)
                    first_tie_index = i
                    tie_flag = 1
                # If no tie, fill in the list
                else:
                    result.append(rank_val)
            # A tie has occured
            else:
                # If another tie occurs
                if curr == next:
                    tie_count = tie_count + 1
                    tie_list.append(rank_val)
                # No more ties, average their values and attach to adjusted list
                else:
                    tie_list.append(rank_val)
                    last_tie_index = i
                    result = _populate_adjusted_vals(tie_list, first_tie_index, \
                        last_tie_index, result)
                    tie_flag = 0
                    tie_count = 0
                    tie_list = []
        
        # If there is a tie that extends to the final position, out of main loop
        # to avoid out of list bounds errors
        if tie_flag == 1:
            tie_list.append(ranks[i+1])
            last_tie_index = i + 1
            result = _populate_adjusted_vals(tie_list, first_tie_index, \
                last_tie_index, result)
        else:
            result.append(ranks[i+1])
            
        return array(result)

    def _populate_adjusted_vals(self, tie_list, first_tie_index, last_tie_index, result):
        """
        Helper function to _remove_ties. Consolidates repeated code
        """
        tie_list = array(tie_list)
        adjusted_value = tie_list.mean()
        while first_tie_index <= last_tie_index:
            result.append(adjusted_value)
            first_tie_index = first_tie_index + 1
        return result

    def _compute_r_value(self, adjusted_ranks, sorted_groups, number_samples):
        """
        Code that performs the actual math involved in solving anosim
        
            PARAMETERS
            adjusted_rank_list: list of the ranks, adjusted for ties
            sorted_group_list: list associating differences to groups
            number_samples: how many total samples
            
            RETURNS
            r: R value computed by anosim
        """
        
        # Compute r_W and r_B
        r_W = mean(adjusted_ranks[sorted_groups==1])
        r_B = mean(adjusted_ranks[sorted_groups==0])
        
        n = number_samples
        divisor = float(n*((n-1)/4))
        r = (r_B - r_W) / divisor 
        return r

    def _format_anosim_results(self, input_path, r_value, p_value='NA'):
        """
        Formats the results of the script for the output file
        
            Line 1: header ("Input filepath ANOSIM_R_VALUE  p_value")
            Line 2: data 
        """
        result = ['Input_filepath\tANOSIM_R_value\tp_value']
        result.append('{0}\t{1}\t{2}'.format(input_path,r_value,p_value))
        return result


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
              file
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


class MantelCorrelogram(CorrelationStats, PermutationStats):
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
                distances between samples (e.g. UniFrac distance matrix)
            geo_dm - a DistanceMatrix object representing some other distance
                measure between samples (most commonly geographical distances,
                but could also be distances in pH, temperature, etc.)
            num_perms - the number of permutations to use when computing the
                p-values
            alpha - the alpha value to use when marking the Mantel
                correlogram plot for significance
        """
        # Can't call super because the two parents take different arguments.
        CorrelationStats.__init__(self, [eco_dm, geo_dm], 2, 3)
        PermutationStats.__init__(self, num_perms)
        self.setAlpha(alpha)

    def getAlpha(self):
        """Returns the alpha value."""
        return self._alpha

    def setAlpha(self, alpha):
        """Sets the alpha value.

        Arguments:
            alpha - the value of alpha. Must be between 0 and 1, inclusive
        """
        if alpha >= 0 and alpha <= 1:
            self._alpha = alpha
        else:
            raise ValueError("Alpha must be between 0 and 1.")

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

                    # Negate the Mantel r statistic because we are using
                    # distance matrices, not similarity matrices (this is a
                    # necessary step, see Legendre's Numerical Ecology
                    # algorithm reference for more details).
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


class Mantel(CorrelationStats, PermutationStats):
    """Class for the Mantel matrix correlation statistical method.

    This class provides the functionality to run a Mantel analysis on two
    distance matrices. A Mantel test essentially computes the Pearson
    correlation between the two distance matrices.
    """

    def __init__(self, dm1, dm2, num_perms, tail_type='two sided'):
        """Constructs a new Mantel instance.

        Arguments:
            dm1 - first DistanceMatrix object to be compared
            dm2 - second DistanceMatrix object to be compared
            num_perms - the number of times to permute the distance matrix
                while calculating the p-value
            tail_type - the type of Mantel test to perform (i.e. hypothesis
                test). Can be "two sided", "less", or "greater"
        """
        CorrelationStats.__init__(self, [dm1, dm2], 2, 3)
        PermutationStats.__init__(self, num_perms)
        self.setTailType(tail_type)

    def runAnalysis(self):
        """Runs a Mantel test over the current distance matrices.

        Returns a dict containing the results. The following keys are set:
            method_name - name of the statistical method
            dm1 - the first DistanceMatrix instance that was used
            dm2 - the second DistanceMatrix instance that was used
            num_perms - the number of permutations used to compute the p-value
            p_value - the p-value computed by the test
            r_value - the Mantel r statistic computed by the test
            perm_stats - a list of Mantel r statistics, one for each
                permutation
            tail_type - the type of Mantel test performed

        Note: R's mantel function will always perform a one-sided test (type
        'greater'), so the p-values may differ from R unless you explicitly
        specify the tail type of 'greater'.
        """
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

        Arguments:
            tail_type - the tail type to use when calculating the p-value.
                Valid types are 'two sided', 'less', or 'greater'.
        """
        if tail_type not in ("two sided", "greater", "less"):
            raise ValueError("Unrecognized alternative hypothesis (tail "
                             "type). Must be either 'two sided', 'greater', "
                             "or 'less'.")
        self._tail_type = tail_type


class PartialMantel(CorrelationStats, PermutationStats):
    """Class for the partial Mantel matrix correlation statistical method.

    This class provides the functionality to run a partial Mantel analysis on
    three distance matrices. A partial Mantel test essentially computes the
    Pearson correlation between two distance matrices after first controlling
    for the effects of a third distance matrix (the control matrix).
    """

    def __init__(self, dm1, dm2, cdm, num_perms):
        """Constructs a new PartialMantel instance.

        Arguments:
            dm1 - first DistanceMatrix object to be compared
            dm2 - second DistanceMatrix object to be compared
            cdm - the control DistanceMatrix object
            num_perms - the number of times to permute while calculating the
                p-value
        """
        CorrelationStats.__init__(self, [dm1, dm2, cdm], 3, 3)
        PermutationStats.__init__(self, num_perms)

    def runAnalysis(self):
        """Runs a partial Mantel test on the current distance matrices.

        Returns a dict containing the results. The following keys are set:
            method_name - name of the statistical method
            mantel_p - the p-value computed by the test
            mantel_r - the Mantel r statistic computed by the test

        Credit: The code herein is based loosely on the implementation found in
        R's vegan package.
        """
        # Calculate the correlation statistic.
        corr = lambda rxy, rxz, ryz: (rxy - rxz*ryz)/(sqrt(1 -
                                      rxz**2)*sqrt(1 - ryz**2))
        # Load initial/placeholder values in the results dictionary.
        res = {}
        res['method_name'] = 'Partial Mantel'
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

        # Calculate the original test statistic (r-value).
        orig_stat = corr(rval1, rval2, rval3)

        # Calculate permuted r-values and p-values, storing
        # them for use in the calculation of the final statistic.
        perm_stats = [0 for i in range(perm_num)]
        numerator = 0
        for i in range(0, perm_num):
            # Permute the first distance matrix and calculate new
            # r and p-values.
            p1 = permute_2d(dm1, permutation(dm1.getSize()))
            dm1_perm = DistanceMatrix(p1, dm1.SampleIds, dm1.SampleIds)
            dm1_perm_flat = dm1_perm.flatten()
            rval1 = pearson(dm1_perm_flat, dm2_flat)
            rval2 = pearson(dm1_perm_flat, cdm_flat)
            perm_stats.append(corr(rval1, rval2, rval3))

            # Sum the permuted statistics for calculation of the final
            # statistic.
            if perm_stats[-1] >= orig_stat:
              numerator += perm_stats[-1]
        # Load the final statistics into the result dictionary.
        res['mantel_r'] = orig_stat
        res['mantel_p'] = (numerator + 1) / (perm_num + 1)
        return res
