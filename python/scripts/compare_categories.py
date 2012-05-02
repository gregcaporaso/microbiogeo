#!/usr/bin/env python
# File created on 22 Apr 2012
from __future__ import division

__author__ = "Michael Dwan"
__copyright__ = "Copyright 2012, The QIIME MiCOS project"
__credits__ = ["Michael Dwan, Logan Knecht"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Michael Dwan"
__email__ = "mdwan.tgen@gmail.com, lgk7@nau.edu"
__status__ = "Development"

from numpy import zeros
from numpy.random import permutation

from os import path, makedirs, listdir

from cogent.util.misc import create_dir

from qiime.format import format_p_value_for_num_iters
from qiime.util import parse_command_line_parameters, make_option
from qiime.parse import parse_distmat, fields_to_dict, \
                        parse_mapping_file, parse_mapping_file_to_dict

#from python.qiime.parse import DistanceMatrix, MetadataMap
#from python.qiime.r_executor import RExecutor
from parse import DistanceMatrix, MetadataMap
from r_executor import RExecutor

#from python.qiime.stats import Anosim
from stats import Anosim, Permanova
#anosim_p_test, _format_anosim_results

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
    # THIS DOESN'T WORK AT ALL AND WON'T CATCH ANY ERRORS FOR A DIR EXISTING
    try:
        if not path.exists(opts.output_dir):
            create_dir(opts.output_dir)
    except:
        option_parser.error("Could not create or access output directory, it already exists. Please delete it and re-run the script"
                                "specified with the -o option.")

    dm_labels, dm_temp = parse_distmat(open(opts.input_dm, 'U'))

    dm = DistanceMatrix(dm_temp, dm_labels, dm_labels)
    md_map = MetadataMap.parseMetadataMap(open(opts.mapping_file))

    cats = opts.categories.split(',')
    first_category = cats[0]

    if   opts.method == 'adonis':
        # verify that category is in mapping file
        map_list = parse_mapping_file(open(opts.mapping_file,'U').readlines())
        if not first_category in map_list[1][1:]:
            print "Category '%s' not found in mapping file columns:" %(first_category)
            print map_list[1][1:]
            exit(1)

        distance_matrix = opts.input_dm
        map_file = opts.mapping_file
        output = opts.output_dir

        command_args = ["-i " + distance_matrix + " -m " + map_file + " -c " + first_category + " -o " + output]

        rex = RExecutor()
        results = rex(command_args, opts.method+".r", output_dir=opts.output_dir)
    elif opts.method == 'anosim':
        anosim_object = Anosim(md_map, dm, first_category)
        runAnalysisOutput = anosim_object(opts.num_permutations)
        outputFile = open(opts.method+"_output_file.txt","w")
        outputFile.write("Method Name:\tR-value:\tP-value:")
        outputFile.write("\n")
        outputFile.write(runAnalysisOutput["method_name"]+"\t"+str(runAnalysisOutput["r_value"])+"\t"+str(runAnalysisOutput["p_value"])+"\t")
        outputFile.write("\n")
        outputFile.close()

    elif opts.method == 'best':
        pass

    elif opts.method == 'dfa':
        #-i otu_table.txt -m map.txt -c BODY_SITE -o dfa_output.txt
        distance_matrix = opts.input_dm
        map_file = opts.mapping_file

        output = opts.output_dir

        command_args = ["-i " + distance_matrix + " -m " + map_file + " -c " + first_category + " -o " + output]

        rex = RExecutor()
        results = rex(command_args, opts.method+".r", output_dir=opts.output_dir)

        print results

        outputFile = open(opts.method+"_output.txt", "w")
        outputFile.write("Method Name:\tR-value:\tP-value:")
        outputFile.write("\n")
        #outputFile.write(runAnalysisOutput["method_name"]+"\t"+str(runAnalysisOutput["r_value"])+"\t"+str(runAnalysisOutput["p_value"])+"\t")
        outputFile.close()
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
        results = rex(command_args, opts.method+".r", output_dir=opts.output_dir)

    elif opts.method == 'mrpp':
        # verify that category is in mapping file
        map_list = parse_mapping_file(open(opts.mapping_file,'U').readlines())
        if not opts.categories in map_list[1][1:]:
            print "Category '%s' not found in mapping file columns:" %(first_category)
            print map_list[1][1:]
            exit(1)

        distance_matrix = opts.input_dm
        map_file = opts.mapping_file
        category = opts.categories
        output = opts.output_dir

        command_args = ["-d " + distance_matrix + " -m " + map_file + " -c " + category + " -o " + output]

        rex = RExecutor()
        results = rex(command_args, opts.method+".r", output_dir=opts.output_dir)
    elif opts.method == 'multicola':
        pass
    elif opts.method == 'permanova':
        #tries to make the file, if it fails it outputs that the file exists and tells the user to delete it before proceeding/
        try:
            outputFile = open(opts.output_dir+"/permanova_output_file.txt", os.O_CREAT)
        except:
            option_parser.error("Could not create the output file. It already exists, delete the file: "+opts.output_dir+"/permanova_output_file.txt in order to resolve this.")
        #makes a permanova object
        permanova_plain = Permanova(md_map, dm, first_category)
        #relies on the __call__ property and
        results = permanova_plain(opts.num_permutations)
        #writes the results to the output dir
        outputFile.write("Method Name:\tR-value:\tP-value:")
        outputFile.write("\n")
        outputFile.write(results["method_name"]+"\t"+str(results["r_value"])+"\t"+str(format_p_value_for_num_iters(results["p_value"], opts.num_permutations))+"\t")
        outputFile.write("\n")
        outputFile.close()
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
        results = rex(command_args, "betadisper.r", output_dir=opts.output_dir)

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
        results = rex(command_args, "rda.r", output_dir=opts.output_dir)

    elif opts.method == 'rm_permanova':
        pass
    else:
        print "Method '%s' not recognized"



if __name__ == "__main__":
    main()
