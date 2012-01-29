#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Andrew Cochran"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Andrew Cochran"]
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "Dan Knights"
__email__ = "danknights@gmail.com"
__status__ = "Development"
 
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from python.qiime.anosim import anosim, anosim_p_test, _format_anosim_results
from qiime.parse import parse_distmat, parse_mapping_file_to_dict
import os

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""ANOSIM"""
script_info['script_description']="""Analysis of similarities (ANOSIM) is a statistical method to test the significance of differences between groups"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""ANOSIM""","""For this script, the user supplies a distance matrix (i.e. resulting file from beta_diversity.py), the output filepath (i.e. results.txt), a mapping file, and the name of the column in the mapping file that contains the grouping information (i.e. "Treatment") as follows:""","""anosim.py -i distmatrix.txt -o result.txt -m mapping.txt -c Treatment"""))
script_info['output_description']="""The resulting output file consists of the R value computed by anosim and the p value (if specified by the original call to the script)"""
script_info['required_options']=[\

make_option('-i', '--input_path',\
     help='path to the input distance matrix file(s) (i.e., the output from beta_diversity.py). [REQUIRED]'),\

make_option('-o', '--output_path',
     help='output path to the name of a single file, [REQUIRED]'),\

make_option('-m', '--map_path',\
     help='path to the location of the mapping file [REQUIRED]'),\

make_option('-c', '--category',\
     help='String which coresponds to the column name containing grouping info [REQUIRED]'),\

]

script_info['optional_options']=[\

make_option('-p', '--ptrials',type='int',\
     help='An integer indicating how many ptrials to be run, [default: 0]',default=0)\

]
script_info['version'] = __version__

def main():
    group_list = {}
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    # Open files
    map_file = open(opts.map_path, 'U')
    input_file = open(opts.input_path,'U')

    # Parse grouping information
    samples, distmtx = parse_distmat(input_file)
    group_maping, comment = parse_mapping_file_to_dict(map_file)
    for sample in group_maping:
        group_list[sample] = group_maping[sample][opts.category]

    # Close files
    map_file.close()
    input_file.close()
    
    outfile = open(opts.output_path, 'w')
    
    # Perform anosim AND p_test
    if opts.ptrials > 0:
        rstat, pvalue  = anosim_p_test(samples, distmtx, group_list, opts.ptrials)
        output = _format_anosim_results(opts.input_path, rstat, pvalue)
        outfile.write('\n'.join(output))
    # Perform ONLY anosim
    else:
        anosim_res = anosim(samples, distmtx, group_list)
        output = _format_anosim_results(opts.input_path, anosim_res)
        outfile.write('\n'.join(output))
        outfile.close()
    map_file.close()

if __name__ == "__main__":
    main()
