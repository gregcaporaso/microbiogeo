#!/usr/bin/env python
from __future__ import division
from numpy import argsort, array, zeros, arange, random, mean, tri, unique
from numpy.random import permutation
import sys
import os.path
from optparse import OptionParser
import cogent.cluster.metric_scaling as ms
from qiime.format import format_coords
from qiime.parse import parse_distmat, parse_mapping_file

__author__ = "Andrew Cochran"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Andrew Cochradifferencesn"]
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "Dan Kgroupingnights"
__email__ = "danknights@gmail.com"
__status__ = "Development"


def anosim(samples, distmtx, group_hash):
    """
    Computes anosim on data provided
    
        PARAMETERS
        samples: list of the sample names
        distmtx: the distance matrix provided
        group_hash: a dict associating the group name to their group
        
        RETURNS
        result: the R value computed by anosim on the data provided 
    """ 
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
    adjusted_rank_list = _remove_ties(sorted_differences, rank_list)
    result = _compute_r_value(adjusted_rank_list, sorted_grouping,\
             matrix_size)
    return result


def anosim_p_test(samples, distmtx, group_hash, ntrials=9999,
        randomfun=random.permutation):
    """
    Performs permutation tests on the provided data
    
        PARAMETERS
        samples: names of the samples
        distmtx: distance matrix provided
        group_hash: each sample associated to a group
        ntrials: how many trials to run
        randomfun: which function to use for randomization (must return value)
        
        RETURNS
        r_value: the R value computed on the data
        p_value: the permutation test result
    """
    # Array to hold value of p-tests
    r_value_permunations = zeros(ntrials)
    
    # Calculate the R value
    r_value = anosim(samples, distmtx, group_hash)

    # Main loop to run the p-tests
    for i in xrange(ntrials):
        
        # Randomize grouping
        grouping_random = []
        for sample in samples:
            grouping_random.append(group_hash[sample])
        grouping_random = randomfun(grouping_random)

        for j, sample in enumerate(samples):
            group_hash[sample] = grouping_random[j]
        r_value_permunations[i] = anosim(samples, distmtx, group_hash)

    # Calculate the p-value and return
    p_value = (sum(r_value_permunations >= r_value) + 1) / (ntrials + 1)
    return r_value, p_value


def _remove_ties(sorted_diffs, ranks):
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


def _populate_adjusted_vals(tie_list, first_tie_index, last_tie_index, result):
    """
    Helper function to _remove_ties. Consolidates repeated code
    """
    tie_list = array(tie_list)
    adjusted_value = tie_list.mean()
    while first_tie_index <= last_tie_index:
        result.append(adjusted_value)
        first_tie_index = first_tie_index + 1
    return result

def _compute_r_value(adjusted_ranks, sorted_groups, number_samples):
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

def _format_anosim_results(input_path, r_value, p_value='NA'):
    """
    Formats the results of the script for the output file
    
        Line 1: header ("Input filepath ANOSIM_R_VALUE  p_value")
        Line 2: data 
    """
    result = ['Input_filepath\tANOSIM_R_value\tp_value']
    result.append('{0}\t{1}\t{2}'.format(input_path,r_value,p_value))
    return result
    
