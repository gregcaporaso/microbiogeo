#!/usr/bin/env python
# File created on 8 Aug 2011
from __future__ import division

__author__ = "Andrew Cochran"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Andrew Cochran"]
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "Dan Knights"
__email__ = "danknights@gmail.com"
__status__ = "Development"

from numpy import max,zeros 
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from microbiogeo.spatial_autocorrelation import morans, morans_variance, _format_morans_results, _w_exponential, _w_inverse
from qiime.parse import parse_distmat, parse_mapping_file_to_dict, parse_mapping_file
import os

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""spatial autocorrelation"""
script_info['script_description']="""spatial autocorrelation is characterized by a multi-dimentional correlation between other nearby locations."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Moran's I""","""For this script, the user supplies a distance matrix (i.e. resulting file from beta_diversity.py), the output filename (i.e. results.txt), a mapping file, and the name of the column in the mapping file that contains vector variable of interest (i.e. "Treatment") as follows:""","""spatial_autocorrelation.py -i distmatrix.txt -o result.txt -m mapping.txt -c Treatment"""))
script_info['output_description']="""The resulting output file consists of the I value computed by spatial autocorrelation and the varience"""
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

make_option('-w', '--weight_method',type='string',\
     help='Which method of weighting the distance matrix values is to be used. Accepts either "inverse" or "exponential" [default: "inverse"]',default="inverse"),\

make_option('-k', '--kvalue',type='int',\
     help='An integer the k value to be used in exponential normalization [default: 1]',default=1),\

make_option('-a', '--variance',type='string',\
     help='Set to true if the varience test should be run, otherwise it will be ignored [default: true]',default="true"),\

make_option('-n', '--normalize_distances',\
     help='If set to true, distances will be normalized to fall between 0 and 1\
     , if set to false they will not be left as is [default: true]',default='true')\


]
script_info['version'] = __version__

def main():
    # Parse command line input
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    # Open files
    map_file = open(opts.map_path, 'U')
    input_file = open(opts.input_path,'U')

    # Parse files
    samples, distmtx = parse_distmat(input_file)
    dict, comment = parse_mapping_file_to_dict(map_file)

    # Extract group column from mapping info
    group_hash = zeros(len(samples))
    for i, sample in enumerate(samples):
        group_hash[i] = dict[sample][opts.category]

    # Close file handles
    map_file.close()
    input_file.close()
    
    # Normalize distances
    if opts.normalize_distances == "true":
        distmtx = distmtx/max(distmtx)
        print "blarg"
        
    # Weight distances
    if opts.weight_method == "inverse":
        w = _w_inverse(distmtx)
    elif opts.weight_method == "exponential":
        w = _w_exponential(distmtx,opts.kvalue)
    else:
        print "Error: Incorrect weight method (-w) specified, use inverse or exponential(must be all lowercase)"
        sys.exit()
    outfile = open(opts.output_path, 'w')
    
    # If varience test is requested
    if opts.variance == "true":
        i_value, v_value  = morans_variance(w, group_hash)
        output = _format_morans_results(opts.input_path, i_value, v_value)
        outfile.write('\n'.join(output))
        
    # Run without varience test
    else:
        morans_res = morans(w, group_hash)
        output = _format_morans_results(opts.input_path, morans_res)
        outfile.write('\n'.join(output))
    outfile.close()
    map_file.close()

if __name__ == "__main__":
    main()
