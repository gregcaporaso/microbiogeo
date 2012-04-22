#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Logan Knecht"
__copyright__ = "Copyright 2012, The QIIME MiCOS Project"
__credits__ = ["Damien Coy, Logan Knecht, Dan Knights"]
__license__ = "GPL"
__version__ = "1.4.0"
__maintainer__ = "Logan Knecht"
__email__ = "lgk7@nau.edu"
__status__ = "Release"
 

from qiime.util import make_option
from os import makedirs, listdir
from os.path import join
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.parse import parse_mapping_file
#from qiime.r_executor import RExecutor
from r_executor import RExecutor
options_lookup = get_options_lookup()
import sys

script_info={}
script_info['brief_description']="""
Analysis of variance using distance matrices - for partitioning distance matrices among sources of variation and fitting linear models to distance matrices.
"""

script_info['script_description']="""
Adonis is a statistical method that operates in the same manner that \
Permutational Multivariate Analysis of Variance(Permanova) does. Adonis is \
 essentially the same thing, except it uses distance matrices. Because it is a \
variation of Permanova, the way which it handles the information passed in is \
similar, except that it can process non-numeric data for analysis as well.

This has been called a "nonparemetric manova". To elaborate on what is being \
said, this variation of the method does not make the same assumptions as MANOVA\
 does. As such it may not determine the same information with the same accuracy.

The primary input that this method is concerned with, is receiving a formula. \
From there the data that is passed in is used in conjunction with the formula \
to calculate the desired information. If desired, you can configure the R \
implementation to accept a formula that uses multiple variables from the \
mapping file passed in.

Outputs: 
    * adonis_results.txt - This file contains the results of the 
          operation. The values outputted to this file are:
        - Df: 
        - SumOfSqs: 
        - MeanSqs: 
        - F.Model: 
        - R2: 
        - Pr(>F)

Source:
http://cran.r-project.org/web/packages/vegan/index.html

This script requires that R be installed and in the search path. To install R \
visit: http://www.r-project.org/. Once R is installed, run R from the \
command-line by typing "R" and pressing enter. Then excecute the command \
"install.packages(vegan)", then type q() to exit."""

script_info['script_usage']=[]

script_info['script_usage'].append(("""Running Adonis on a distance matrice and its related map file. This uses the PH column from map.txt as the category.""",""" """,""" python morans_i.py -i unweighted_unifrac_dm.txt -m map.txt -c PH -o moransIOutput """))

script_info['output_description']="""Outputs a file that contains the values \
used to calculate the p-value for Adonis.
"""

script_info['required_options'] = [\
    make_option('-i', '--input_data', help='This is the argument used for the \
distance matrix being used. It should correspond to the mapping file being \
passed in.'),
    make_option('-m', '--mapping_file', help='This is the corresponding mapping\
 file of the distance matrix being used. This is where the category is selected\
 from.'),
    make_option('-c', '--category', help='Name of meta data category selected \
from the mapping file in order to identify the spatial correlation based on \
that category.'),
]

script_info['optional_options']=[\
    make_option('-o','--output_dir',default='.',\
            help='the output directory [default: %default]'),
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

    command_args = ["-i " + distance_matrix + " -m " + map_file + " -c " + category + " -o " + output]

    rex = RExecutor()
    results = rex(command_args, "adonis.r", output_dir=opts.output_dir, remove_tmp=True)
        
if __name__ == "__main__":
    main()
