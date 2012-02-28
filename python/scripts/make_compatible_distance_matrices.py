#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jr378@nau.edu"
__status__ = "Development"
 
from qiime.format import format_distance_matrix
from qiime.parse import parse_distmat
from qiime.util import make_compatible_distance_matrices, make_option, \
    parse_command_line_parameters

script_info = {}
script_info['brief_description'] = """
Makes two distance matrices compatible based on sample IDs
"""
script_info['script_description'] = """
This script finds the intersection of two distance matrices and creates two \
resulting distance matrices containing that intersection, with sample IDs in \
the same order. An example usage of this script is prepping two distance \
matrices for use in a Mantel-type test that takes two distance matrices and \
computes their correlation. Many of the functions that do this in R assume \
that each corresponding cell in the two matrices are for the same pairwise \
comparison (i.e. cell (1, 10) in dm 1 contains distances between sample 1 and \
samle 10, so cell (1, 10) in dm 2 is assumed to contain the distance between \
sample 1 and sample 10).

Please note that this functionality is already used in the existing Mantel \
test in QIIME (compare_distance_matrices.py), but this script is provided for \
convenience when using an external tool, such as R, to do the processing.
"""
script_info['script_usage'] = [("Make two distance matrices compatible",
    "This example shows how to create two new distance matrices that are "
    "compatible.",
    "%prog -i unweighted_unifrac_dm.txt,PH_dm.txt -o "
    "unweighted_unifrac_dm_comp.txt,PH_dm_comp.txt")]
script_info['output_description'] = """
The output is two distance matrices that only contain sample IDs that were in \
both input distance matrices, and the sample IDs are in the same order for \
both output matrices.
"""
script_info['required_options'] = [
    make_option('-i', '--input_dms',
        help='the two input distance matrices, comma-separated'),
    make_option('-o', '--output_dms',
        help='filepaths for the two output distance matrices, '
             'comma-separated')]
script_info['optional_options'] = []
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # Open the input distance matrices, parse them, find the intersection, and
    # write the two new distance matrices to the output filepaths.
    input_dm_fps = opts.input_dms.split(',')
    output_dm_fps = opts.output_dms.split(',')
    if len(input_dm_fps) != 2 or len(output_dm_fps) != 2:
        option_parser.error("You must provide exactly two input and output "
            "distance matrix filepaths.")

    labels1, dm1_data = parse_distmat(open(input_dm_fps[0], 'U'))
    labels2, dm2_data = parse_distmat(open(input_dm_fps[1], 'U'))

    (dm1_labels, dm1), (dm2_labels, dm2) = make_compatible_distance_matrices(
        parse_distmat(open(input_dm_fps[0],'U')),
        parse_distmat(open(input_dm_fps[1],'U')))
    assert (dm1_labels == dm2_labels), "The order of sample IDs is not the " +\
        "same for the two matrices."

    output1_f = open(output_dm_fps[0], 'w')
    output2_f = open(output_dm_fps[1], 'w')
    output1_f.write(format_distance_matrix(dm1_labels, dm1))
    output2_f.write(format_distance_matrix(dm2_labels, dm2))
    output1_f.close()
    output2_f.close()

if __name__ == "__main__":
    main()
