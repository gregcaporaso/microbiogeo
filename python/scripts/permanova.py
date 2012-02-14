#!/usr/bin/env python
# File created on 18 Jul 2011
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
from python.qiime.permanova import permanova, permanova_p_test,\
     _format_permanova_results
from qiime.parse import parse_distmat, parse_mapping_file_to_dict,\
parse_mapping_file
import os

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""permanova"""
script_info['script_description']="""PERMANOVA: permutational multivariate analysis of variance on the basis of any distance measure"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""PERMANOVA""","""For this script, the user supplies a distance matrix (i.e. resulting file from beta_diversity.py), the output filename (i.e. results.txt), a mapping file, and the name of the column in the mapping file that contains the grouping information (i.e. "Treatment") as follows:""","""permanova.py -i distmatrix.txt -o result.txt -m mapping.txt -c Treatment"""))
script_info['output_description']="""The resulting output file consists of the F value computed by permanova and the p value (if specified by the original call to the script)"""
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
    # Local Vars
    group_hash = {}
    
    # Parse command line input
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    # Open files
    map_file = open(opts.map_path, 'U')
    input_file = open(opts.input_path,'U')

    # Parse files
    samples, distmtx = parse_distmat(input_file)
    group_mapping, comment = parse_mapping_file_to_dict(map_file)
    
    # Extract group column from mapping info
    for sample in group_mapping:
        group_hash[sample] = group_mapping[sample][opts.category]

    # Close file handles
    map_file.close()
    input_file.close()
    
    outfile = open(opts.output_path, 'w')
    
    # If permunation test is requested
    if opts.ptrials > 0:
        rstat, pvalue  = permanova_p_test(samples, distmtx, group_hash, opts.ptrials)
        output = _format_permanova_results(opts.input_path, rstat, pvalue)
        outfile.write('\n'.join(output))
    # Run without permunation test
    else:
        permanova_res = permanova(samples, distmtx, group_hash)
        output = _format_permanova_results(opts.input_path, permanova_res)
        outfile.write('\n'.join(output))
        outfile.close()
    map_file.close()

if __name__ == "__main__":
    main()
