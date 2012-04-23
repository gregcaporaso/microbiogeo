#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Damien Coy"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Damien Coy"]
__license__ = "GPL"
__version__ = "1.4.0"
__maintainer__ = "Damien Coy"
__email__ = "damien.coy@nau.edu"
__status__ = "Release"

from qiime.util import make_option
from os import makedirs, listdir
from os.path import join
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.parse import parse_mapping_file
from qiime.r_executor import RExecutor
options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Run PERMDISP (betadisper) on given distance matrix using map and category."""
script_info['script_description']="""This script utilizes the r_executor class to call the/
 r implementation of PERMDISP (betadisper) using the specified distance matrix, otu table, category, and/
 output file provided. The results are then written to the specified file."""

script_info['script_usage']=[]
script_info['script_usage'].append(("""Simple example""","""""","""permdisp.py -i distance_matrix.txt -m map.txt -c Category -o result"""))
script_info['output_description']="""Outputs the results of the PERMDISP (betadisper)"""
script_info['required_options'] = [\
    make_option('-i', '--input_data', help='The input data file containing the distance matrix'),
    make_option('-m', '--mapping_file', help='The mapping file that corresponds to the distance matrix'),
    make_option('-c', '--category', help='The name of the category in the mapping file'),
    make_option('-o', '--output_dir', help='The output directory that will contain the results'),
]
script_info['optional_options']=[\
    make_option('-f','--force',action='store_true',\
        dest='force',help='Force overwrite of existing output directory'+\
        ' (note: existing files in output_dir will not be removed)'+\
        ' [default: %default]'),
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # create the output directories
    try:
        makedirs(opts.output_dir)
    except OSError:
        if opts.force:
            pass
        else:
            # This check helps users avoid overwriting previous output.
            print "Output directory already exists. Please choose "+\
             "a different directory, or force overwrite with -f."
            exit(1)

    # verify that category is in mapping file
    map_list = parse_mapping_file(open(opts.mapping_file,'U').readlines())
    if not opts.category in map_list[1][1:]:
        print "Category '%s' not found in mapping file columns:" %(opts.category)
        print map_list[1][1:]
        exit(1)

    distance_matrix = opts.input_data
    map_file = opts.mapping_file
    category = opts.category
    output = opts.output_dir

    command_args = ["-d " + distance_matrix + " -m " + map_file + " -c " + category + " -o " + output]

    rex = RExecutor()
    results = rex(command_args, "betadisper.r", output_dir=opts.output_dir, remove_tmp=True)

if __name__ == "__main__":
    main()