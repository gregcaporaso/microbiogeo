#!/usr/bin/env python
# File created on 22 Apr 2012
from __future__ import division

__author__ = "Michael Dwan"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Michael Dwan"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Michael Dwan"
__email__ = "mdwan.tgen@gmail.com"
__status__ = "Development"


from os import path, makedirs, listdir

from cogent.util.misc import create_dir

from qiime.util import parse_command_line_parameters, make_option
from qiime.parse import parse_distmat, fields_to_dict, \
                        parse_mapping_file

from python.qiime.parse import DistanceMatrix, MetadataMap
from python.qiime.r_executor import RExecutor

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [\
 # All methods use these
make_option('--method', help='The category analysis method. Valid options: [adonis, anosim, bioenv, dfa, isa, lsa, morans_i, mrpp, multicola, permanova, permdisp, rda, rm_permanova]'),\
 make_option('-i','--input_dm',help='the input distance matrix'),\
 make_option('-o','--output_dir',help='the output directory [default: %default]', default='.'),\
 make_option('-m','--mapping_file', help='Mapping file'),
 make_option('-c','--categories',help='A comma delimited list of categories from the mapping file(NOTE: many methods take just a single category, if multiple are passed only the first will be selected.)'),\
]
script_info['optional_options'] = [\
 # All methods use these
 make_option('-n','--num_permutations',help='the number of iterations to perform',default=100,type='int'),
]
script_info['version'] = __version__



def main():
    option_parser, opts, args =\
    parse_command_line_parameters(**script_info)

    # Create the output dir if it doesn't already exist.
    try:
        if not path.exists(opts.output_dir):
            create_dir(opts.output_dir)
    except:
        option_parser.error("Could not create or access output directory "
                                "specified with the -o option.")


    dm_labels, dm_temp = parse_distmat(open(opts.input_dm, 'U'))

    dm = DistanceMatrix(dm_temp, dm_labels, dm_labels)
    md_map = MetadataMap.parseMetadataMap(open(opts.mapping_file))

    cats = opts.categories.split(',')

    if   opts.method == 'adonis':
        pass
    elif opts.method == 'anosim':
        pass
    elif opts.method == 'best':
        pass
    elif opts.method == 'dfa':
        pass
    elif opts.method == 'isa':
        pass
    elif opts.method == 'lsa':
        pass
    elif opts.method == 'morans_i':
        category = cats[0]
        # verify that category is in mapping file
        map_list = parse_mapping_file(open(opts.mapping_file,'U').readlines())
        if not category in map_list[1][1:]:
            print "Category '%s' not found in mapping file columns:" %(category)
            print map_list[1][1:]
            exit(1)

        distance_matrix = opts.input_dm
        map_file = opts.mapping_file
        category = category
        output = opts.output_dir

        command_args = ["-i " + distance_matrix + " -m " + map_file + " -c " + category + " -o " + output]

        rex = RExecutor()
        results = rex(command_args, "morans_i.r", output_dir=opts.output_dir, remove_tmp=True)

    elif opts.method == 'mrpp':
        pass
    elif opts.method == 'multicola':
        pass
    elif opts.method == 'permanova':
        pass

    elif opts.method == 'permdisp':
        category = cats[0]
        # verify that category is in mapping file
        map_list = parse_mapping_file(open(opts.mapping_file,'U').readlines())
        if not category in map_list[1][1:]:
            print "Category '%s' not found in mapping file columns:" %(category)
            print map_list[1][1:]
            exit(1)

        distance_matrix = opts.input_dm
        map_file = opts.mapping_file
        output = opts.output_dir

        command_args = ["-i " + distance_matrix + " -m " + map_file + " -c " + category + " -o " + output]

        rex = RExecutor()
        results = rex(command_args, "betadisper.r", output_dir=opts.output_dir, remove_tmp=True)

    elif opts.method == 'rda':
        category = cats[0]
        # verify that category is in mapping file
        map_list = parse_mapping_file(open(opts.mapping_file,'U').readlines())
        if not category in map_list[1][1:]:
            print "Category '%s' not found in mapping file columns:" %(category)
            print map_list[1][1:]
            exit(1)

        distance_matrix = opts.input_dm
        map_file = opts.mapping_file
        output = opts.output_dir

        command_args = ["-i " + distance_matrix + " -m " + map_file + " -c " + category + " -o " + output]

        rex = RExecutor()
        results = rex(command_args, "rda.r", output_dir=opts.output_dir, remove_tmp=True)

    elif opts.method == 'rm_permanova':
        pass
    else:
        print "Method '%s' not recognized"


    # 'adonis'
    # 'anosim'
    # 'best'
    # 'dfa'
    # 'isa'
    # 'lsa'
    # 'morans_i'
    # 'mrpp'
    # 'multicola'
    # 'permanova'
    # 'permdisp'
    # 'rda'
    # 'rm_permanova'


if __name__ == "__main__":
    main()