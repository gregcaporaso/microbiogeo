#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Logan Knecht"
__copyright__ = "Copyright 2011, The QIIME MiCOS Project"
__credits__ = ["Damien Coy, Logan Knecht"]
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
This function computes Moran.s I autocorrelation coefficient of x giving a \
matrix of weights using the method described by Gittleman and Kot (1990).
"""
script_info['script_description']="""
This script performs a Moran's I calculation on a distance matrix and its \
mapping file. The Moran's I calculation identifies the spatial correlation \
of a set of data based on the category selected. This category must be a \
column listed in the mapping file. Using that selection it then calculates \
the information between the mapping file and the distance matrix. Then it \
returns an output value between negative one (-1) to positive one (1). When the\
 result is closer to a negative one (-1) that means that the population is \
dispersed. When the result is closer to a positive (1) it means that there is \
clustering of those samples spatially. When the output is closer to zero (0) \
that indicates that there is no discernable pattern being observed from the \
data.

Outputs: 
    * Morans_I_results.txt - This file contains the results of the Moran's I
          operation. There are four values outputted to this file.
        - observed: Moran's I index of x.
        - expected: Expected value of I under the null hypothesis.
        - sd: The standard deviation of I under the null hypothesis.
        - p.value: The p-value of having the observed value under the null \
hypothesis.

Source:
http://cran.r-project.org/web/packages/ape/ape.pdf

Example:
python morans_i.py -i unweighted_unifrac_dm.txt -m map.txt -c PH -o moransIOutput

This script requires that R be installed and in the search path. To install R \
visit: http://www.r-project.org/. Once R is installed, run R from the \
command-line by typing "R" and pressing enter. Then excecute the command \
"install.packages(ape)", then type q() to exit."""

script_info['script_usage']=[]

script_info['script_usage'].append(("""Running a Moran's I calculation using PH as the selected category from map.txt""",""" """,""" python morans_i.py -i unweighted_unifrac_dm.txt -m map.txt -c PH -o moransIOutput """))

script_info['output_description']="""Outputs a ranking of features (e.g. OTUs) by importance, an estimation of the generalization error of the classifier, and the predicted class labels and posterior class probabilities \
according to the classifier."""

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

#Not sure what this is about - LK
#errortype_choices = ['oob','loo','cv5','cv10']

script_info['optional_options']=[\
    make_option('-o','--output_dir',default='.',\
            help='the output directory [default: %default]'),
    make_option('-f','--force',action='store_true',\
        dest='force',help='Force overwrite of existing output directory'+\
        ' (note: existing files in output_dir will not be removed)'+\
        ' [default: %default]'),
    make_option('--ntree',type='int',default=500,\
        help='Number of trees in forest (more is better but slower) [default: %default]'),
    #make_option('-e', '--errortype',type='choice',default='oob',
        #choices = errortype_choices,
        #help='type of error estimation. Valid choices are: ' +\
            #', '.join(errortype_choices) + '. '+\
            #'oob: out-of-bag, fastest, only builds one classifier, use for quick estimates; ' +\
            #'cv5: 5-fold cross validation, provides mean and standard deviation of error, use for good estimates on very large data sets; ' +\
            #'cv10: 10-fold cross validation, provides mean and standard deviation of error, use for best estimates; ' +\
            #'loo: leave-one-out cross validation, use for small data sets (less than ~30-50 samples) ' +\
            #'[default %default]')
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
    results = rex(command_args, "morans_i.r", output_dir=opts.output_dir, remove_tmp=True)
        
if __name__ == "__main__":
    main()
