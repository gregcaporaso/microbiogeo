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
 
from random import shuffle
from qiime.format import format_distance_matrix
from qiime.parse import parse_distmat
from qiime.util import parse_command_line_parameters, make_option

script_info = {}
script_info['brief_description'] = "Shuffles the labels of a distance matrix"
script_info['script_description'] = """
This script shuffles the sample ID labels of a distance matrix. This \
functionality may be useful, for example, if you have a distance matrix that \
you are using a positive control in testing and you want to create a negative \
control to use in additional testing.
"""
script_info['script_usage'] = [("Shuffle distance matrix labels",
    "This example shows how to shuffle the labels of a distance matrix.",
    "%prog -i unweighted_unifrac_dm.txt -o "
    "shuffled_unweighted_unifrac_dm.txt")]
script_info['output_description'] = """
The output is a distance matrix with the same data as the input distance \
matrix, but with the labels shuffled.
"""
script_info['required_options'] = [
    make_option('-i','--input_distance_matrix',
        help='the input distance matrix', type="existing_filepath"),
    make_option('-o', '--output_distance_matrix',
        help='path to store the output distance matrix', type="new_filepath")]
script_info['optional_options'] = []
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # Open the input distance matrix and parse it. Shuffle its labels and write
    # them and the original data to the output file.
    labels, dm_data = parse_distmat(open(opts.input_distance_matrix, 'U'))
    shuffle(labels)
    output_f = open(opts.output_distance_matrix, 'w')
    output_f.write(format_distance_matrix(labels, dm_data))
    output_f.close()

if __name__ == "__main__":
    main()
