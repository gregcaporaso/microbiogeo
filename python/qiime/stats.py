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

from python.qiime.parse import DistanceMatrix, MetadataMap


class DistanceMatrixStats(object):
    """Base class for distance matrix-based statistical methods.

    This class provides an interface to setting and accessing an arbitrary
    number of distance matrices. Users of this class can optionally specify the
    number of allowable distance matrices and their minimum allowable size (the
    default is no restrictions on either of these).

    It is the parent class of CorrelationStats and CategoryStats.
    """

    def __init__(self, dms, num_dms=-1, min_dm_size=-1):
        """Default constructor.

        Initializes an instance with the provided list of distance matrices.

        Arguments:
            dms - a list of DistanceMatrix objects
            num_dms - the exact number of allowable distance matrices. If -1
                (the default), there is no restriction on how many distance
                matrices the user can set
            min_dm_size - the minimum size that all distance matrices must have
                that are stored by this instance. If -1, no size restriction
        """
        self._num_dms = num_dms
        self._min_dm_size = min_dm_size
        self.setDistanceMatrices(dms)

    def getDistanceMatrices(self):
        """Returns the list of distance matrices."""
        return self._dms

    def setDistanceMatrices(self, dms):
        """Sets the list of distance matrices to the supplied list.

        Arguments:
            dms - the new list of distance matrices being assigned
        """
        if not isinstance(dms, ListType):
            raise TypeError("The item passed in as the new list was not a "
                            "list data type.")
        if self._num_dms >= 0 and len(dms) != self._num_dms:
            raise ValueError("Cannot set %d distance matrices. Must provide "
                             "exactly %d distance matrices." % (len(dms),
                             self._num_dms))
        for dm in dms:
            if not isinstance(dm, DistanceMatrix):
                raise TypeError('Invalid type: %s; expected DistanceMatrix' %
                                dm.__class__.__name__)
            elif self._min_dm_size >= 0 and dm.getSize() < self._min_dm_size:
                raise ValueError("Distance matrix of size %dx%d is smaller "
                                 "than the minimum allowable distance matrix "
                                 "size of %dx%d for this analysis." %
                                 (dm.getSize(), dm.getSize(),
                                  self._min_dm_size, self._min_dm_size))
        self._dms = dms

    def __call__(self, num_perms=999):
        """Runs the statistical method and returns relevant results.

        The return value of this method is a python dictionary with arbitrary
        key/value pairs of results, since each statistical method returns
        different results.

        This method returns an empty result set (it is essentially not
        implemented) and should be implemented by subclasses to perform their
        specific statistical analysis. Subclasses should call the parent
        class' __call__ method first to obtain any results from the parent and
        then add more results to the dict that is obtained from the parent.

        Arguments:
            num_perms - the number of permutations to use in the statistical
                method. If the method is not permutation-based, simply ignore
                this argument
        """
        if num_perms < 0:
            raise ValueError("The number of permutations must be greater than "
                             "or equal to zero.")
        return {}


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

    def __init__(self, dms, num_dms=-1, min_dm_size=-1):
        """Default constructor.

        Creates a new instance with the provided list of distance matrices.

        Arguments:
            dms - a list of DistanceMatrix objects
            num_dms - the exact number of allowable distance matrices. If -1
                (the default), there is no restriction on how many distance
                matrices the user can set
            min_dm_size - the minimum size that all distance matrices must have
                that are stored by this instance. If -1, no size restriction
        """
        super(CorrelationStats, self).__init__(dms, num_dms, min_dm_size)

    def setDistanceMatrices(self, dms):
        """Sets the list of distance matrices to the supplied list.

        This method overrides the parent method and enforces more checks to
        ensure that at least one distance matrix is provided and that all of
        the distance matrices are compatible.

        Arguments:
            dms - the new list of distance matrices being assigned
        """
        super(CorrelationStats, self).setDistanceMatrices(dms)
        if len(dms) < 1:
            raise ValueError("Must provide at least one distance matrix.")

        size = dms[0].getSize()
        sample_ids = dms[0].SampleIds
        for dm in dms:
            if dm.getSize() != size:
                raise ValueError("All distance matrices must have the same "
                                 "number of rows and columns.")
            if dm.SampleIds != sample_ids:
                raise ValueError("All distance matrices must have matching "
                                 "sample IDs.")


class CategoryStats(DistanceMatrixStats):
    """Base class for categorical statistical analyses.

    It is subclassed by categorical statistical methods such as DB-RDA or BEST.
    Categorical statistical methods usually have some categorical grouping of
    samples, and the significance of this grouping is usually what is tested.
    For example, are treatment samples significantly different from control
    samples? This is not always the case (e.g. DB-RDA is an ordination
    technique), but most of the categorical methods follow this general design.

    A valid instance of CategoryStats must have at least one distance matrix
    and a single metadata map containing the sample IDs of the distance matrix
    or matrices.
    """

    def __init__(self, mdmap, dms, cats, num_dms=-1, min_dm_size=-1):
        """Default constructor.

        Creates a new instance with the provided distance matrices,
        metadata map, and list of categories.

        Arguments:
            mdmap - a MetadataMap instance
            dms - a list of DistanceMatrix objects
            cats - a list of strings denoting categories in the metadata map
                that will be used by this analysis (i.e. the grouping
                variable(s))
            num_dms - the exact number of allowable distance matrices. If -1
                (the default), there is no restriction on how many distance
                matrices the user can set
            min_dm_size - the minimum size that all distance matrices must have
                that are stored by this instance. If -1, no size restriction
        """
        super(CategoryStats, self).__init__(dms, num_dms, min_dm_size)
        self.setMetadataMap(mdmap)
        self.setCategories(cats)
        self._validate_compatibility()

    def setMetadataMap(self, new_mdmap):
        """Sets the instance's metadata map.

        Arguments:
            new_mdmap - A MetadataMap object instance
        """
        if not isinstance(new_mdmap, MetadataMap):
            raise TypeError('Invalid type: %s; not MetadataMap' %
                            new_mdmap.__class__.__name__)
        self._metadata_map = new_mdmap

    def getMetadataMap(self):
        """Returns the instance's metadata map.

        The metadata map is returned as a MetadataMap class instance.
        """
        return self._metadata_map

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

    def _validate_compatibility(self):
        """Raises an error if the current dms and map are incompatible.

        This method will raise an error if any of the sample IDs in any of the
        distance matrices are not found in the metadata map. Ordering of
        sample IDs is not taken into account.
        """
        for dm in self.getDistanceMatrices():
            for samp_id in dm.getSampleIds():
                if samp_id not in self.getMetadataMap().getSampleIds():
                    raise ValueError("The sample ID '%s' was not found in the "
                                     "metadata map." % samp_id)

    def __call__(self, num_perms=999):
        """Runs the statistical method and returns relevant results.

        The return value of this method is a python dictionary with arbitrary
        key/value pairs of results, since each statistical method returns
        different results.

        This method returns an empty result set (it is essentially not
        implemented) and should be implemented by subclasses to perform their
        specific statistical analysis. Subclasses should call the parent
        class' __call__ method first to obtain any results from the parent and
        then add more results to the dict that is obtained from the parent.

        Arguments:
            num_perms - the number of permutations to use in the statistical
                method. If the method is not permutation-based, simply ignore
                this argument
        """
        # Make sure the *current* distance matrices and metadata map are
        # compatible before continuing.
        self._validate_compatibility()
        return super(CategoryStats, self).__call__(num_perms)


class Anosim(CategoryStats):
    """Class for the ANOSIM categorical statistical analysis.

    Briefly, ANOSIM tests whether two or more groups of samples are
    significantly different. The user of the class specifies a category in the
    metadata map to group samples by.

    This code is heavily based on Andrew Cochran's original procedural version.
    """

    def __init__(self, mdmap, dm, cat, random_fn=permutation):
        """Initializes an instance with the specified analysis parameters.

        Arguments:
            mdmap - the MetadataMap instance to obtain grouping info from
            dm - the DistanceMatrix instance to obtain distances from
            cat - the category string to group samples by (must be in the
                metadata map)
            random_fn - the function to use when randomizing the grouping
                during calculation of the p-value. It must return a value and
                must be callable
        """
        super(Anosim, self).__init__(mdmap, [dm], [cat], num_dms=1)
        self.setRandomFunction(random_fn)

    def getRandomFunction(self):
        """Returns the randomization function used in p-value calculations."""
        return self._random_fn

    def setRandomFunction(self, random_fn):
        """Setter for the randomization function used in p-value calcs.

        Arguments:
            random_fn - the function to use when randomizing the grouping
                during calculation of the p-value. It must return a value and
                must be callable
        """
        if hasattr(random_fn, '__call__'):
            self._random_fn = random_fn
        else:
            raise TypeError("The supplied function reference is not callable.")

    def __call__(self, num_perms=999):
        """Runs ANOSIM on the current distance matrix and sample grouping.

        Returns a dict containing the results. The following keys are set:
            method_name - name of the statistical method
            r_value - the ANOSIM R statistic computed by the test
            p_value - the p-value computed by the test, or 'NA' if the number
                of permutations was zero

        Arguments:
            num_perms - the number of permutations to use when calculating the
                p-value
        """
        results = super(Anosim, self).__call__(num_perms)
        category = self.getCategories()[0]
        samples = self.getDistanceMatrices()[0].getSampleIds()

        # Create the group map, which maps sample ID to category value (e.g.
        # sample 1 to 'control' and sample 2 to 'fast').
        group_map = {}
        for samp_id in samples:
            group_map[samp_id] = self.getMetadataMap().getCategoryValue(
                    samp_id, category)

        # Calculate the R statistic with the grouping found in the current
        # metadata map.
        r_stat = self._anosim(group_map)

        if num_perms > 0:
            # Calculate the p-value based on the number of permutations.
            perm_stats = []
            for i in range(num_perms):
                # Randomize grouping. We don't use values() in order to
                # preserve ordering in case the user's random function doesn't
                # change the order of the items in the list.
                grouping_random = [group_map[sample] for sample in samples]
                grouping_random = self.getRandomFunction()(grouping_random)
                for j, sample in enumerate(samples):
                    group_map[sample] = grouping_random[j]
                perm_stats.append(self._anosim(group_map))
            # Calculate the p-value.
            p_value = (sum(perm_stats >= r_stat) + 1) / (num_perms + 1)
        else:
            p_value = 'NA'

        results['method_name'] = 'ANOSIM'
        results['r_value'] = r_stat
        results['p_value'] = p_value
        return results

    def _anosim(self, group_map):
        """Computes ANOSIM on the supplied grouping, returning the R value.

        The R value is between -1 and 1 and indicates the strength of the
        grouping.

        Arguments:
            group_map - a python dict mapping sample ID to category value (e.g.
                sample 1 to 'control' and sample 2 to 'fast'). This map must
                contain a key for each sample ID in the current distance
                matrix
        """
        dm = self.getDistanceMatrices()[0]
        dm_size = dm.getSize()

        # Create grouping matrix, where a one means that the two samples are in
        # the same group (e.g. control) and a zero means that they aren't.
        within_between = zeros((dm_size, dm_size))
        for i, i_sample in enumerate(dm.getSampleIds()):
            for j, j_sample in enumerate(dm.getSampleIds()):
                if group_map[i_sample] == group_map[j_sample]:
                    within_between[i][j] = 1

        # Extract upper triangle from the distance and grouping matrices.
        distances = dm.getDataMatrix()[tri(dm_size) == 0]
        grouping = within_between[tri(dm_size) == 0]

        # Sort extracted data.
        sorted_distances = []
        sorted_grouping = []
        for idx in argsort(distances):
             sorted_distances.append(distances[idx])
             sorted_grouping.append(grouping[idx])

        # Account for rank ties, then compute R statistic.
        rank_list = range(1, len(sorted_distances) + 1)
        adjusted_rank_list = self._remove_ties(sorted_distances, rank_list)
        return self._compute_r_value(adjusted_rank_list, sorted_grouping,
                                     dm_size)

    def _remove_ties(self, sorted_dists, ranks):
        """Replaces repeat values with the average of them.

        Returns a list containing the adjusted ranks.

        Arguments:
            sorted_dists: list of the sorted distances
            ranks: list containing the ranks of each of the differences
        """
        result = []
        ties = []
        tie_count = 0
        tie_flag = 0

        for i in range(len(sorted_dists) - 1):
            # Store state information.
            curr_dist = sorted_dists[i]
            next_dist = sorted_dists[i+1]
            rank_val = ranks[i]

            # A tie has not occured yet.
            if tie_flag == 0:
                if curr_dist == next_dist:
                    # We have a tie, so add the current rank to the tie list.
                    tie_count = tie_count + 1
                    ties.append(rank_val)
                    first_tie_index = i
                    tie_flag = 1
                else:
                    # If no tie, fill in the list with the current rank.
                    result.append(rank_val)
            else:
                # A tie has already occured.
                if curr_dist == next_dist:
                    # If another tie occurs, add the current rank to the tie
                    # list.
                    tie_count = tie_count + 1
                    ties.append(rank_val)
                else:
                    # No more ties, average their values and attach to adjusted
                    # rank list.
                    ties.append(rank_val)
                    last_tie_index = i
                    result.extend(self._get_adjusted_vals(ties,
                            first_tie_index, last_tie_index))
                    tie_flag = 0
                    tie_count = 0
                    ties = []
        # If there is a tie that extends to the final position, we must process
        # it here to avoid out of list bounds errors.
        if tie_flag == 1:
            ties.append(ranks[i+1])
            last_tie_index = i + 1
            result.extend(self._get_adjusted_vals(ties, first_tie_index,
                                                  last_tie_index))
        else:
            result.append(ranks[i+1])
        return result

    def _get_adjusted_vals(self, ties, first_tie_idx, last_tie_idx):
        """Helper function to _remove_ties. Consolidates repeated code."""
        adjusted_val = sum(ties) / len(ties)
        return [adjusted_val] * ((last_tie_idx - first_tie_idx) + 1)

    def _compute_r_value(self, adjusted_ranks, sorted_groups, num_samps):
        """Code that performs the actual math involved in solving ANOSIM.

        Returns the ANOSIM R value (between -1 and 1).

        Arguments:
            adjusted_ranks - list of the ranks, adjusted for ties
            sorted_groups - list associating distances to groups
            num_samps: how many total samples
        """
        adjusted_ranks = array(adjusted_ranks)
        sorted_groups = array(sorted_groups)

        # Compute r_W and r_B.
        r_W = mean(adjusted_ranks[sorted_groups==1])
        r_B = mean(adjusted_ranks[sorted_groups==0])
        divisor = num_samps * ((num_samps - 1) / 4)
        return (r_B - r_W) / divisor

class Permanova(CategoryStats):
    """This code is heavily based on Andrew Cochran's original procedural version."""

    def __init__(self, mdmap, dm, cat, random_fn=permutation):
        """Initializes an instance with the specified analysis parameters.

        Arguments:
            mdmap - the MetadataMap instance to obtain grouping info from
            dm - the DistanceMatrix instance to obtain distances from
            cat - the category string to group samples by (must be in the
                metadata map)
            num_perms - the number of permutations to use when calculating the
                p-value. If zero, the p-value will not be calculated. Must be
                greater than or equal to zero
            random_fn - the function to use when randomizing the grouping
                during calculation of the p-value. It must return a value and
                must be callable
        """
        super(Permanova, self).__init__(mdmap, [dm], [cat], num_dms=1)
        self.setRandomFunction(random_fn)

    def getRandomFunction(self):
        """Returns the randomization function used in p-value calculations."""
        return self._random_fn

    def setRandomFunction(self, random_fn):
        """Setter for the randomization function used in p-value calcs.

        Arguments:
            random_fn - the function to use when randomizing the grouping
                during calculation of the p-value. It must return a value and
                must be callable
        """
        if hasattr(random_fn, '__call__'):
            self._random_fn = random_fn
        else:
            raise TypeError("The supplied function reference is not callable.")

    def _permanova(self, samples, distmtx, grouping):
        """Computes PERMANOVA pseudo-f-statistic

           PARAMETERS
       grouping: a Metamap object
        """


        # Local Vars
        unique_n = []       # number of samples in each group
        group_map = {}
        map = {}
        metaMap = self.getMetadataMap()

        #make map
        for sample in metaMap.getSampleIds():
                subkey = {}
                for cat in metaMap.getCategoryNames():
                    subkey[cat] = metaMap.getCategoryValue(sample, cat)
                map[sample] = subkey


        # Extract the unique list of group labels
        gl_unique = unique(array(grouping.values()))

        # Calculate number of gorups and unique 'n's
        number_groups = len(gl_unique)
        for i, i_string in enumerate(gl_unique):
            group_map[i_string] = i
            unique_n.append(grouping.values().count(i_string))

        # Create grouping matrix
        grouping_matrix = -1 * ones((len(distmtx),len(distmtx)))
        for i, i_sample in enumerate(samples):
            grouping_i = grouping[i_sample]
            for j, j_sample in enumerate(samples):
                if grouping_i == grouping[j_sample]:
                    grouping_matrix[i][j] = group_map[grouping[i_sample]]

        # Extract upper triangle
        distances = distmtx[tri(len(distmtx)) == 0]
        gropuing = grouping_matrix[tri(len(grouping_matrix)) == 0]

        # Compute f value
        result = self._compute_f_value(distances,gropuing,len(distmtx),number_groups,unique_n)
        return result


    def permanova_p_test(self, samples, distmtx, group_list, ntrials=9999,\
                     randomfun=random.permutation):
        """Performs the calculations for the permutation test

        PARAMETERS
        samples: names of the samples
        distmtx: the data
        group_list: listing of the grouping information
        ntrials: how many trials to run, default 9999

        RETURNS
        f_value: the value of permanova
        p_value: permutation factor

        """
        # Array to store permutation values
        f_value_permunations = zeros(ntrials)

        # Calculate the F-Value
        f_value = self._permanova(samples,distmtx,group_list)

        # Run p-tests
        for i in xrange(ntrials):

            # Randomize Grouping
            grouping_random = []
            for sample in samples:
                grouping_random.append(group_list[sample])
            grouping_random = randomfun(grouping_random)

            # Calculate p-values
            for j, sample in enumerate(samples):
                group_list[sample] = grouping_random[j]
            f_value_permunations[i] = self._permanova(samples,distmtx,group_list)
       
        p_value = (sum(f_value_permunations >= f_value) + 1) / (ntrials + 1)
        return f_value, p_value


    def _compute_f_value(self, distances, groupings, number_samples, number_groups, unique_n):
        """Performs the calculations for the f value

           PARAMETERS
           difference_list: a list of the distance values
           group_list: a list associating the distances to their groups
           number_samples: how many samples there are
           number_groups: how many groups there are
           unique_n: list containing how many samples are in each within group
        """
        a = number_groups                 # number of groups
        N = number_samples                # total samples

        # Calculate s_T
        s_T = sum(distances*distances)/N

        # Calculate s_W for each group, this accounts for diff group sizes
        s_W = 0
        for i in range(number_groups):
            group_ix = groupings==i
            diffs = distances[group_ix]
            s_W = s_W + sum(diffs**2)/unique_n[i]

        # Execute the formula
        s_A = s_T - s_W
        f = (s_A/(a-1))/(s_W/(N-a))
        return f


class BioEnv(CategoryStats):
    """Class for BioEnv analysis."""

    def __init__(self, dm, metadata_map, cats):
        """Default constructor."""

        super(BioEnv, self).__init__(metadata_map, [dm], cats, num_dms=1)

    def __call__(self, num_perms=999):
        """Runs the BioEnv analysis on a distance matrix using specified
        metadata map categories.

        TODO: ADD COMMENTS

        Returns a dictionary which contains the resulting data. Keys:
            method_name - name of the statistical method
        """

        res = super(BioEnv, self).__call__(num_perms)
        cats = self.getCategories()
        dm = self.getDistanceMatrices()[0]
        dm_flat = dm.flatten()
        dm_flat_ranked = self._get_rank(dm_flat)

        row_count = dm.getSize()
        col_count = len(cats)
        sum = 0
        stats = [(-777, '') for c in range(col_count+1)]
        for i in range(col_count+1):
            if i < 11:
                combo = list(combinate([j for j in range(0,col_count)], i))[1:]
            else:
                combo = list(combinate([j for j in range(0,col_count+1)], i))[0:1]

            for c in range(len(combo)):
                cat_mat = self._make_cat_mat(cats, combo[c])
                cat_dm = self._derive_euclidean_dm(cat_mat, row_count)
                cat_dm_flat_ranked = self._get_rank(cat_dm.flatten())
                r = self._spearman_correlation(dm_flat_ranked, cat_dm_flat_ranked, ranked=True)
                if r > stats[i][0]:
                    stats[i] = (r, ','.join(str(s+1) for s in combo[c]))

        # for s in stats[1:]:
        #     print s
        # print (2**col_count - 1)/2

        res['method_name'] = 'BioEnv'
        res['num_vars'] = col_count
        res['vars'] = ['%s = %d' % (name,val+1) for val,name in enumerate(cats)]
        res['bioenv_rho_vals'] = stats[1:]

        return res


    def _derive_euclidean_dm(self, cat_mat, dim):
        """Returns an n x n, euclidean distance matrix, where n = len(cats) """

        dm_labels = self.getDistanceMatrices()[0].getSampleIds()
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

        dm = self.getDistanceMatrices()[0]
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

#if __name__ == '__main__':
    # dm = DistanceMatrix.parseDistanceMatrix(open('unweighted_unifrac_dm.txt'))
    # md_map = MetadataMap.parseMetadataMap(open('vars.txt'))
#    dm = DistanceMatrix.parseDistanceMatrix(open('dm.txt'))
#    md_map = MetadataMap.parseMetadataMap(open('vars2.txt'))
#
#
#    cats = ['TOT_ORG_CARB', 'SILT_CLAY', 'ELEVATION', 'SOIL_MOISTURE_DEFICIT', 'CARB_NITRO_RATIO', 'ANNUAL_SEASON_TEMP', 'ANNUAL_SEASON_PRECPT', 'PH', 'CMIN_RATE', 'LONGITUDE', 'LATITUDE']
#
#    bioenv = BioEnv(dm, md_map, cats)
    # a = [1,2,4,3,1,6,7,8,10,4]
    # b = [2,10,20,1,3,7,5,11,6,13]
    # x = (1,  2, 4, 3, 1, 6, 7, 8, 10, 4, 100, 2, 3, 77)
    # y = (2, 10, 20, 1, 3, 7, 5, 11, 6, 13, 5, 6, 99, 101)
    # r = (1,2,4,5,2,2,4,3,1,4)
    # s = (2,3,5,4,2,2,3,4,3,2)
    # u = (1,2,3,4,5,6,7,8,9)
    # v = (10, 11, 4, 2, 9, 33, 1, 5, 88)
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

        mdmap = self.getMetadataMap()
        dm = self.getDistanceMatrices()[0]
        k = dm.getSize() - 1
        inertia_str = "Distance squared"

        if dm.getMax() >= 4:
            inertia_str = "mean " + inertia_str
            adjust = 1
        else:
            adjust = sqrt(k)

        points, eigs = principal_coordinates_analysis(dm.getDataMatrix())

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

    def __init__(self, eco_dm, geo_dm, alpha=0.05):
        """Constructs a new MantelCorrelogram instance.

        Arguments:
            eco_dm - a DistanceMatrix object representing the ecological
                distances between samples (e.g. UniFrac distance matrix)
            geo_dm - a DistanceMatrix object representing some other distance
                measure between samples (most commonly geographical distances,
                but could also be distances in pH, temperature, etc.)
            alpha - the alpha value to use when marking the Mantel
                correlogram plot for significance
        """
        super(MantelCorrelogram, self).__init__([eco_dm, geo_dm], num_dms=2,
                                                min_dm_size=3)
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

    def __call__(self, num_perms=999):
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

        Arguments:
            num_perms - the number of permutations to use when calculating the
                p-values

        Note: This code is heavily based on the implementation of
        mantel.correlog in R's vegan package.
        """
        results = super(MantelCorrelogram, self).__call__(num_perms)
        eco_dm = self.getDistanceMatrices()[0]
        geo_dm = self.getDistanceMatrices()[1]
        dm_size = eco_dm.getSize()

        # Find the number of lower triangular elements (excluding the
        # diagonal).
        num_dists = dm_size * (dm_size - 1) // 2

        # Use Sturge's rule to determine the number of distance classes.
        num_classes = int(ceil(1 + log2(num_dists)))

        # Create the matrix of distance classes. Each element in the matrix
        # contains what distance class the original element is in. Also find
        # the distance class indices, which are the midpoints in each distance
        # class.
        dist_class_matrix, class_indices = self._find_distance_classes(geo_dm,
            num_classes)

        # Start assembling the results.
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
                                         tail_type='greater')
                    mantel_test_results = mantel_test(num_perms)
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
                        p_val = perm_sum / (num_perms + 1)
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
    distance matrices. A Mantel test essentially computes the Pearson
    correlation between the two distance matrices.
    """

    def __init__(self, dm1, dm2, tail_type='two sided'):
        """Constructs a new Mantel instance.

        Arguments:
            dm1 - first DistanceMatrix object to be compared
            dm2 - second DistanceMatrix object to be compared
            tail_type - the type of Mantel test to perform (i.e. hypothesis
                test). Can be "two sided", "less", or "greater"
        """
        super(Mantel, self).__init__([dm1, dm2], num_dms=2, min_dm_size=3)
        self.setTailType(tail_type)

    def __call__(self, num_perms=999):
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

        Arguments:
            num_perms - the number of times to permute the distance matrix
                while calculating the p-value

        Note: R's mantel function will always perform a one-sided test (type
        'greater'), so the p-values may differ from R unless you explicitly
        specify the tail type of 'greater'.
        """
        resultsDict = super(Mantel, self).__call__(num_perms)
        results = self._mantel_test(num_perms)

        resultsDict['method_name'] = "Mantel"
        resultsDict['dm1'] = self.getDistanceMatrices()[0]
        resultsDict['dm2'] = self.getDistanceMatrices()[1]
        resultsDict['num_perms'] = num_perms
        resultsDict['p_value'] = results[0]
        resultsDict['r_value'] = results[1]
        resultsDict['perm_stats'] = results[2]
        resultsDict['tail_type'] = self.getTailType()

        return resultsDict

    def _mantel_test(self, n):
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

        Arguments:
            n - the number of permutations
        """
        m1, m2 = self.getDistanceMatrices()
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


class PartialMantel(CorrelationStats):
    """Class for the partial Mantel matrix correlation statistical method.

    This class provides the functionality to run a partial Mantel analysis on
    three distance matrices. A partial Mantel test essentially computes the
    Pearson correlation between two distance matrices after first controlling
    for the effects of a third distance matrix (the control matrix).
    """

    def __init__(self, dm1, dm2, cdm):
        """Constructs a new PartialMantel instance.

        Arguments:
            dm1 - first DistanceMatrix object to be compared
            dm2 - second DistanceMatrix object to be compared
            cdm - the control DistanceMatrix object
        """
        super(PartialMantel, self).__init__([dm1, dm2, cdm], num_dms=3,
                                            min_dm_size=3)

    def __call__(self, num_perms=999):
        """Runs a partial Mantel test on the current distance matrices.

        Returns a dict containing the results. The following keys are set:
            method_name - name of the statistical method
            mantel_p - the p-value computed by the test
            mantel_r - the Mantel r statistic computed by the test

        Arguments:
            num_perms - the number of times to permute the distance matrix
                while calculating the p-value

        Credit: The code herein is based loosely on the implementation found in
        R's vegan package.
        """
        res = super(PartialMantel, self).__call__(num_perms)

        # Calculate the correlation statistic.
        corr = lambda rxy, rxz, ryz: (rxy - rxz*ryz)/(sqrt(1 -
                                      rxz**2)*sqrt(1 - ryz**2))
        # Load initial/placeholder values in the results dictionary.
        res['method_name'] = 'Partial Mantel'
        res['mantel_r'] = None
        res['mantel_p'] = None

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
        perm_stats = [0 for i in range(num_perms)]
        numerator = 0
        for i in range(0, num_perms):
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
        res['mantel_p'] = (numerator + 1) / (num_perms + 1)
        return res
