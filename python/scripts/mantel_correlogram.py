#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Rideout"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Jai Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

from os import path
from cogent.util.misc import create_dir
from qiime.format import format_p_value_for_num_iters
from qiime.parse import parse_distmat, fields_to_dict
from qiime.util import (get_options_lookup, make_compatible_distance_matrices,
                        make_option, parse_command_line_parameters)
from python.qiime.parse import DistanceMatrix
from python.qiime.stats import MantelCorrelogram

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """
Computes a Mantel correlogram between two distance matrices
"""
script_info['script_description'] = """
This script creates a Mantel correlogram, which is a plot of distance classes \
versus Mantel statistics. Briefly, an ecological distance matrix (e.g. \
UniFrac distance matrix) and a second distance matrix (e.g. spatial \
distances, pH distances, etc.) are provided to the script. The second \
distance matrix has its distances split into a number of distance classes \
(the number of classes is determined by Sturge's rule). A Mantel test is run \
over these distance classes versus the ecological distance matrix. The \
Mantel statistics obtained from each of these tests is then plotted in a \
correlogram. A filled-in point on the plot indicates that the Mantel \
statistic was statistically significant (you may provide what alpha to use).
"""
script_info['script_usage'] = [("Compute Mantel correlogram",
"This example computes a Mantel correlogram on two distance matrices using "
"999 permutations in each Mantel test. Output is written to the "
"mantel_output directory.",
"%prog -i unweighted_unifrac_dm.txt,weighted_unifrac_dm.txt -o mantel_output "
"-n 999")]
script_info['output_description']= """
Two files are created in the output directory: a text file containing \
information about the distance classes, their associated Mantel statistics \
and p-values, etc. and an image of the correlogram plot.
"""
script_info['required_options'] = [
    make_option('-i', '--input_dms',
        help='the two input distance matrices, comma-separated'),
    options_lookup['output_dir']
]
script_info['optional_options'] = [
    make_option('-n', '--num_permutations',
        help='the number of permutations to perform', default=100, type='int'),
    make_option('-a', '--alpha',
        help='the value of alpha to use when denoting significance in the '
             'correlogram plot', default=0.05, type='float'),
    make_option('-g', '--image_type',
        help='type of image to produce (i.e. png, svg, pdf) '
             '[default: %default]', default='pdf', type="choice",
        choices=['pdf', 'png', 'svg']),
    make_option('-s', '--sample_id_map_fp',
        help='Map of original sample ids to new sample ids '
             '[default: %default]', default=None)
]
script_info['version'] = __version__

comment = """# Number of entries refers to the number of rows (or cols) 
# retained in each distance matrix after filtering the distance matrices 
# to include only those samples that were in both distance matrices. 
# p-value contains the correct number of significant digits.
# Distance classes with values of None were in the second half of the distance
# classes and not all samples could be included in the distance class, so
# calculations were not performed.
"""

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    num_perms = opts.num_permutations
    alpha = opts.alpha

    # Create the output dir if it doesn't already exist.
    try:
        create_dir(opts.output_dir)
    except:
        option_parser.error("Could not create or access output directory "
                            "specified with the -o option.")
    input_dm_fps = opts.input_dms.split(',')
    results_f = open(path.join(
        opts.output_dir, "mantel_correlogram_results.txt"), 'w')
    results_f.write(comment)

    # Make the two distance matrices compatible before running the analysis.
    # This code was taken from Greg's compare_distance_matrices.py script.
    sample_id_map_fp = opts.sample_id_map_fp
    if sample_id_map_fp:
        sample_id_map = dict([(k,v[0]) \
         for k,v in fields_to_dict(open(sample_id_map_fp, "U")).items()])
    else:
        sample_id_map = None

    (dm1_labels, dm1), (dm2_labels, dm2) = make_compatible_distance_matrices(
        parse_distmat(open(input_dm_fps[0], 'U')),
        parse_distmat(open(input_dm_fps[1], 'U')), lookup=sample_id_map)
    if len(dm1_labels) < 3:
        option_parser.error("The distance matrices were not large enough "
            "after filtering them to include only samples that were in both "
            "matrices. The minimum required size to compute a Mantel "
            "correlogram is 3x3.")



    # Write header info to the results file.
    results_f.write('DM1: %s\nDM2: %s\nNumber of entries: %d\nNumber of '
        'permutations: %d\nAlpha: %s\n' % (input_dm_fps[0], input_dm_fps[1],
        len(dm1_labels), num_perms, alpha))

    # Construct a MantelCorrelogram object and run the analysis.
    results = MantelCorrelogram(DistanceMatrix(dm1, dm1_labels,
        dm1_labels), DistanceMatrix(dm2, dm2_labels, dm2_labels),
        alpha=alpha)(num_perms)

    # Write the correlogram plot to a file.
    results['correlogram_plot'].savefig(path.join(opts.output_dir,
        'mantel_correlogram.%s' % opts.image_type), format=opts.image_type)
    
    # Iterate over the results and write them to the text file.
    results_f.write('\nClass index\tNum dists\tMantel stat\tp-val\t'
        'p-val (Bonferroni corrected)\n')
    for class_idx, num_dist, r, p, p_corr in zip(results['class_index'],
        results['num_dist'], results['mantel_r'], results['mantel_p'],
        results['mantel_p_corr']):
        if p is not None:
            p_str = format_p_value_for_num_iters(p, num_perms)
        else:
            p_str = None
        if p_corr is not None:
            p_corr_str = format_p_value_for_num_iters(p_corr, num_perms)
        else:
            p_corr_str = None
        results_f.write('%s\t%d\t%s\t%s\t%s\n' % (class_idx, num_dist, r, p,
            p_corr))
    results_f.close()


if __name__ == "__main__":
    main()
