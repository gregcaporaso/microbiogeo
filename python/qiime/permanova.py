#!/usr/bin/env python
# File created on 18 Jul 2011

from __future__ import division
from numpy import tri, array, ones, random, unique, zeros
from numpy.random import permutation
import sys
import os.path
from optparse import OptionParser
import cogent.cluster.metric_scaling as ms
from qiime.format import format_coords
from qiime.parse import parse_distmat, parse_mapping_file

__author__ = "Andrew Cochran"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Andrew Cochran"]
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "Dan Knights"
__email__ = "danknights@gmail.com"
__status__ = "Development"


def permanova(samples, distmtx, grouping):
    """Computes PERMANOVA pseudo-f-statistic
    
       PARAMETERS
       samples: list of the sample labels
       distmtx: a distance matrix as returned by parse_distmat
       grouping: a dict, keys = sample labels, values = which group they're in
    """
    # Local Vars
    group_map = {}      # dict to map group number to group name
    unique_n = []       # number of samples in each group
    
    # Extract the unique list of group labels 
    gl_unique = unique(array(grouping.values()))
 
    # Calculate number of gorups and unique 'n's
    number_groups = len(gl_unique)
    for i, i_string in enumerate(gl_unique):
        group_map[i_string] = i
        unique_n.append(grouping.values().count(i_string))

    # Create grouping matrix
    grouping_matrix = -1 * ones(distmtx.shape)
    for i, i_sample in enumerate(samples):
        grouping_i = grouping[i_sample]
        for j, j_sample in enumerate(samples):
            if grouping_i == grouping[j_sample]:
                grouping_matrix[i][j] = group_map[grouping[i_sample]]

    # Extract upper triangle
    distances = distmtx[tri(len(distmtx)) == 0]
    gropuing = grouping_matrix[tri(len(grouping_matrix)) == 0]

    # Compute f value
    result = _compute_f_value(distances,gropuing,len(distmtx),number_groups,unique_n)
    return result


def permanova_p_test(samples, distmtx, group_list, ntrials=9999,\
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
    f_value = permanova(samples, distmtx, group_list)

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
        f_value_permunations[i] = permanova(samples, distmtx, group_list)

    p_value = (sum(f_value_permunations >= f_value) + 1) / (ntrials + 1)
    return f_value, p_value


def _compute_f_value(distances, groupings, number_samples, number_groups, unique_n):
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

def _format_permanova_results(input_path, r_value, p_value='NA'):
    """Formats the data to be output to a text file
    """
    result = ['Input_filepath\tpermanova_R_value\tp_value']
    result.append('{0}\t{1}\t{2}'.format(input_path,r_value,p_value))
    return result
 
