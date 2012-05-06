#!/usr/bin/env python
# File created on 22 Apr 2012
"""
This is a file that aggregates the compare category methods and provides \
a consistent interface for using these statistical methods. 
"""
from __future__ import division

__author__ = "Logan Knecht"
__copyright__ = "Copyright 2012, The QIIME MiCOS project"
__credits__ = ["Logan Knecht, Michael Dwan"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Logan Knecht"
__email__ = "lgk7@nau.edu, mdwan.tgen@gmail.com"
__status__ = "Development"

from os import path, makedirs, listdir

from numpy import zeros
from numpy.random import permutation

from cogent.util.misc import create_dir

from qiime.format import format_p_value_for_num_iters
from qiime.util import parse_command_line_parameters, make_option
from qiime.parse import parse_distmat, fields_to_dict, \
                        parse_mapping_file, parse_mapping_file_to_dict

#from python.qiime.parse import DistanceMatrix, MetadataMap
#from python.qiime.r_executor import RExecutor
from parse import DistanceMatrix, MetadataMap
from qiime.r_executor import RExecutor

#from python.qiime.stats import Anosim, Permanova
from stats import Anosim, BioEnv, Permanova

script_info = {}
script_info['brief_description'] = """
Analyzes distance matrices for information using statistical methods
"""
script_info['script_description'] = """
This script allows for the anaylsis of distance matrices, (in the case of the \
dfa method it's an otu table), using several statistical methods. These methods\
 are Adonis, Anosim, BEST, DFA, ISA(not implemented), LSA(not implemented), \
Moran's I, MRPP, Multicola(not implemented), PERMANOVA, PERMDISP, RDA, and \
Repeated Measures PERMANOVA(not implemented).

Adonis - This method takes a distance matrix and mapping file. It then \
identifies important points in the data and performs F-tests on the initial \
data, and random permutations(via  shuffling) the category data. Then, \
it finally returns the information that was identified in the samples. It's \
stated that it partitions (or seperates the data) for this analysis in order\
 to find underlying relationship.

Anosim - This method takes in a distance matrix and a mapping file. \
This method tests whether two or more categories are significantly \
different. You can specify a category in the metadata mapping file to separate \
samples into groups and then test whether there are significant differences \
between those groups. For example, you might test whether Control samples are \
significantly different from Fast samples. Since ANOSIM is non-parametric, \
significance is determined through permutations.

BEST - This method looks at the numerical environmental variables relating \
samples in a distance matrix. For instance, the unifrac distance matrix and \
pH and latitude (or any other number of variables) in soil samples, and ranks \
them in order of which best explain patterns in the communities.

DFA - This method takes in an OTU table and a mapping file. This method \
is heavily related to ANOVA/MANOVA, except that it operates in the opposite \
direction. ANOVA/MANOVA test whether categorical independent \
variable(s) effectively predict continuous dependent variable(s). On the other \
hand, DFA tests how well one or more independent continuous variables predicts \
a categorical dependent variable. In a nutshell, DFA seeks to answer the \
question 'Is this set of variables effective at predicting category/group \
membership?'. When applied to ecology, DFA can answer the question 'Are these \
species/OTUs effective at predicting category/group membership?' or 'What is \
the best combination of species/OTUs for predicting category/group membersip?'

Moran's I - This method takes in a distance matrix and mapping file. Then it \
uses the geographical data to identify what type of spatial configuration the \
samples have. Are they dispersed, clustered, or of no distinctly noticeable \
configuration when compared to each other?

MRPP - This method takes in a distance matrix and a mapping file. It then \
tests whether two or more categories are significantly different. You can \
specify a category in the metadata mapping file to separate samples into \
groups and then test whether there are significant differences between those \
groups. For example, you might test whether Control samples are significantly \
different from Fast samples. Since MRPP is non-parametric, significance is \
determined through permutations.

PERMANOVA - This method takes distance matrix and a mapping file. This method \
is for testing the simultaneous response of one or more variables to one or \
more factors in an ANOVA experimental design on the basis of any distance \
metric. It returns a R value and a P value. The first thing it does is \
calculate the distances between each pair of sampled units to obtain a \
distance matrix. It then calculates the test-statistics from this according \
to the relevant experimental design.

PERMDISP - This method takes a distance matrix and a mapping file. \
This method is a procedure for the analysis of multivariate homogeneity of \
group dispersions (variances). Permutations can be utilized to measure the \
dissimilatities between groups.

RDA - This method takes a distance matrix and a mapping file. This method \
is an ordination method that shows grouping/clustering of samples based on \
a category in the metadata mapping file and a distance matrix. This category \
is used to explain the variability between samples. Thus, RDA is similar to \
PCoA except that it is constrained, while PCoA is unconstrained (you must \
specify which category should be used to explain the variability in your data).
"""

script_info['script_usage'] = []

script_info['script_usage'].append(("Adonis",
"Performs the Adonis statistical method on a distance matrix and mapping file "
"using the HOST_SUBJECT_ID category and 999 permutations. Then it outputs the "
"results to the 'adonis' directory. The full file path will be: "
"./adonis/adonis_results.txt",
"%prog --method adonis -i datasets/keyboard/unweighted_unifrac_dm.txt -m \
datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o adonis -n 999"))

script_info['script_usage'].append(("Anosim",
"Performs the Anosim statistical method on a distance matrix and mapping file "
"using the HOST_SUBJECT_ID category and 999 perutations. Then it outputs the \
results to the 'anosim' directory. The full file path will be: \
./anosim/anosim_results.txt",
"%prog --method anosim -i datasets/keyboard/unweighted_unifrac_dm.txt -m \
datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o anosim -n 999"))

script_info['script_usage'].append(("BEST",
"Performs the BEST statistical method on a distance matrix and mapping file "
"using the LATITUDE and LONGITUDE categories. Then it outputs the results to \
the 'best' directory. The full file path will be: ./best/best_results.txt",
"%prog --method best -i datasets/keyboard/unweighted_unifrac_dm.txt -m \
datasets/keyboard/map.txt -c LATITUDE,LONGITUDE -o best"))

#-------------------------------------------------------------------------
#TODO FIX THIS LATER SOMETHING SOMETHING BIOM FILE MESSING WITH THIS
script_info['script_usage'].append(("DFA",
"Performs the DFA statistical method on a distance matrix and mapping file "
"using the HOST_SUBJECT_ID category. Then it outputs the results to \
the 'dfa' directory. The full file path will be: ./dfa/dfa_results.txt",
"%prog --method dfa  -i datasets/keyboard/unweighted_unifrac_dm.txt -m \
datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o dfa"))
#-------------------------------------------------------------------------

script_info['script_usage'].append(("Moran's I",
"Performs the Moran's I statistical method on a distance matrix and mapping \
file using the PH category. Then it outputs the results to the 'morans_i' \
directory. The full file path will be: ./morans_i/Morans_I_results.txt",
"%prog --method morans_i -i  datasets/88_soils/unweighted_unifrac_dm.txt -m \
datasets/88_soils/map.txt -c PH -o morans_i"))

script_info['script_usage'].append(("MRPP", "Performs the MRPP statistical \
method on a distance matrix and mapping file using the HOST_SUBJECT_ID \
category. Then it outputs the results to the 'mrpp' directory. The full file \
path will be: ./mrpp/mrpp_results.txt",
"%prog --method mrpp -i datasets/keyboard/unweighted_unifrac_dm.txt -m \
datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o mrpp -n 999"))

script_info['script_usage'].append(("Multicola", "", ""))

script_info['script_usage'].append(("PERMANOVA", "Performs the PERMANOVA \
statistical method on a distance matrix and mapping file using the \
HOST_SUBJECT_ID category. Then it outputs the results to the 'permanova' \
directory. The full file path will be: ./permanova/permanova_results.txt",
"%prog --method permanova -i datasets/keyboard/unweighted_unifrac_dm.txt -m \
datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o permanova -n 999"))

script_info['script_usage'].append(("PERMDISP", "Performs the PERMDISP \
statistical method on a distance matrix and mapping file using the \
HOST_SUBJECT_ID category. Then it outputs the results to the 'permdisp' \
directory. The full file path will be: ./permdisp/betadisper_results.txt",
"%prog --method permdisp -i datasets/keyboard/unweighted_unifrac_dm.txt -m \
datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o permdisp -n 999"))

script_info['output_description']= """
Adonis: 
One file is created and outputs the results into it. The results will be:\
Analysis of variance(AOV) table, degrees of freedom, sequential sums of \
squares, mean squares, F statistics, partial R-squared and p-values, based \
on the N permutations.

Anosim: 

Best: 

DFA: 

Moran's I: 

MRPP: 

Multicola: 

PERMANOVA: 

PERMDISP: 

RDA: 

"""

script_info['required_options'] = [\
 # All methods use these
make_option('--method', help='The category analysis method. Valid options: \
    [adonis, anosim, best, dfa, morans_i, mrpp, multicola, \
    permanova, permdisp, rda]'),\
 make_option('-i','--input_dm',help='This should be a distance matrix that is \
being passed in, unless the method being performed is DFA. If that is the case \
the DFA method requires that an otu table be passed in instead.'),\
 make_option('-o','--output_dir',help='the output directory \
 [default: %default]', default='.'),\
 make_option('-m','--mapping_file', help='Mapping file'),
 make_option('-c','--categories',help='A comma delimited list of categories \
     from the mapping file(NOTE: many methods take just a single category, if\
     multiple are passed only the first will be selected.)'),\
]
script_info['optional_options'] = [\
 # All methods use these
 make_option('-n','--num_permutations',help='the number of iterations to \
     perform',default=999,type='int'),
]
script_info['version'] = __version__

def main():
    """
    This is the entry point for the script to run.
    """
    option_parser, opts, args =\
    parse_command_line_parameters(**script_info)

    # Create the output dir if it doesn't already exist.
    #TODO THIS DOESN'T WORK AT ALL AND WON'T CATCH ANY ERRORS FOR A DIR EXISTING
    try:
        if not path.exists(opts.output_dir):
            create_dir(opts.output_dir)
    except:
        option_parser.error("Could not create or access output directory, it \
            already exists. Please delete it and re-run the script"
                                "specified with the -o option.")

    #parses dist mat for all methods not DFA
    if opts.method != 'dfa':
        dm_labels, dm_temp = parse_distmat(open(opts.input_dm, 'U'))
        dm = DistanceMatrix(dm_temp, dm_labels, dm_labels)

    #parse mapping file
    md_map = MetadataMap.parseMetadataMap(open(opts.mapping_file))

    #separates all categerios into a list, then grabs the first category
    categories = opts.categories.split(',')
    first_category = categories[0]
    
    #cursory check to make sure all categories passed in are in mapping file
    maps = parse_mapping_file(open(opts.mapping_file,'U').readlines())
    for category in categories:
        if not category in maps[1][1:]:
            print "Category '%s' not found in mapping file columns:" % category
            print maps[1][1:]
            exit(1)

    if opts.method == 'adonis':
        command_args = ["-d " + opts.input_dm + " -m " + opts.mapping_file + \
            " -c " + first_category + " -o " + opts.output_dir + " -n " + \
            str(opts.num_permutations)]
        rex = RExecutor()
        rex(command_args, "adonis.r", output_dir=opts.output_dir)
    elif opts.method == 'anosim':
        #runs anosim
        anosim = Anosim(md_map, dm, first_category)
        anosim_results = anosim(opts.num_permutations)
        #anosim has been run, now writing results to file
        output_file = open(opts.output_dir + "/" + opts.method + \
            "_results.xt","w+")
        output_file.write("Method Name:\tR-value:\tP-value:")
        output_file.write("\n")
        output_file.write(anosim_results["method_name"]+"\t"+\
            str(anosim_results["r_value"])+"\t"+\
            str(anosim_results["p_value"])+"\t")
        output_file.write("\n")
        output_file.close()
    elif opts.method == 'best':
        #makes a bioenv object 
        bioenv = BioEnv(dm, md_map, categories)
        #relies on the __call__ property and returns the results
        bioenv_results = bioenv()
        #writes the output to a file
        output_file = open(opts.output_dir+"/best_results.txt", 'w+')
        output_file.write("Method Name:\tNum_Vars:\t")
        output_file.write("\n")
        output_file.write(bioenv_results["method_name"]+"\t"+\
            str(bioenv_results["num_vars"]) + "\t")
        output_file.write("\n")
        output_file.write("Variables:\t")
        output_file.write("\n")
        for variable in bioenv_results["vars"]:
            output_file.write(str(variable) + "\t")
        output_file.write("\n")
        output_file.write("RHO_Values:\t")
        output_file.write("\n")
        for rho_val in bioenv_results["bioenv_rho_vals"]:
            output_file.write(str(rho_val) + "\t")
        output_file.write("\n")
        output_file.close()
    elif opts.method == 'dfa':
        #-----------------------------------------------------------------
        #TODO TRY TO SEE IF THIS WILL WORK FOR MULTIPLE CATEGORIES
        command_args = ["-i " + opts.input_dm + " -m " + opts.mapping_file + \
            " -c " + first_category + " -o " + opts.output_dir]
        #command_args = ["-i " + opts.input_dm + " -m " + opts.mapping_file + \
            #" -c " + "".join(categories) + " -o " + opts.output_dir]
        rex = RExecutor()
        rex(command_args, "dfa.r", output_dir=opts.output_dir)
        #-----------------------------------------------------------------
    elif opts.method == 'morans_i':
        command_args = ["-i " + opts.input_dm + " -m " + opts.mapping_file + \
            " -c " + first_category + " -o " + opts.output_dir]
        rex = RExecutor()
        rex(command_args, "morans_i.r", output_dir=opts.output_dir)
    elif opts.method == 'mrpp':
        command_args = ["-d " + opts.input_dm + " -m " + opts.mapping_file + \
            " -c " + first_category + " -o " + opts.output_dir + \
            " -n " + str(opts.num_permutations)]
        rex = RExecutor()
        rex(command_args, "mrpp.r", output_dir=opts.output_dir)
    elif opts.method == 'multicola':
        pass
    elif opts.method == 'permanova':
        #makes a permanova object
        permanova_plain = Permanova(md_map, dm, first_category)
        #relies on the __call__ property and returns the results
        permanova_results = permanova_plain(opts.num_permutations)
        #writes the results to the output dir
        output_file = open(opts.output_dir+"/permanova_results.txt", 'w+')
        output_file.write("Method Name:\tR-value:\tP-value:")
        output_file.write("\n")
        output_file.write(permanova_results["method_name"]+"\t"+\
            str(permanova_results["r_value"]) + "\t" + \
            format_p_value_for_num_iters(permanova_results["p_value"], \
            opts.num_permutations)+"\t")
        output_file.write("\n")
        output_file.close()
    elif opts.method == 'permdisp':
        command_args = ["-d " + opts.input_dm + " -m " + opts.mapping_file + \
            " -c " + first_category + " -o " + opts.output_dir]
        rex = RExecutor()
        rex(command_args, "betadisper.r", output_dir=opts.output_dir)
    elif opts.method == 'rda':
        command_args = ["-i " + opts.input_dm + " -m " + opts.mapping_file + \
            " -c " + first_category + " -o " + opts.output_dir]
        rex = RExecutor()
        rex(command_args, "rda.r", output_dir=opts.output_dir)
    else:
        print "Method '%s' not recognized" % opts.method

if __name__ == "__main__":
    main()
