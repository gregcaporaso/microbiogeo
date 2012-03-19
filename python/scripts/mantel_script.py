import sys
sys.path.append("/home/ubuntu/")
sys.path.append("/home/ubuntu/biom-format-0.9.1/biom-format-0.9.1/python-code/biom/")

from qiime.util import make_option
from qiime.parse import parse_distmat
from qiime.format import format_p_value_for_num_iters
from qiime.util import (parse_command_line_parameters, 
                        get_options_lookup,
                        make_compatible_distance_matrices)

from numpy import array, asarray, ravel, sqrt
from numpy.random import permutation

from python.qiime.mantel import Mantel

import re

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso, Logan Knecht"]
__license__ = "GPL"
__version__ = "1.4.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Script for computing Mantel correlations between as set of distance matrices"
script_info['script_description'] = ""
script_info['script_usage'] = [("","Perform Mantel test on all pairs of four distance matrices, including 1000 Monte Carlo iterations. Write the output to mantel_out.txt.","%prog -i weighted_unifrac_dm.txt,unweighted_unifrac_dm.txt,weighted_unifrac_even100_dm.txt,unweighted_unifrac_even100_dm.txt -o mantel_out.txt -n 1000")]
script_info['output_description']= ""

script_info['required_options'] = [\
# Example required option
make_option('-i','--input_dms',help='the input distance matrices, comma-separated'),\
make_option('-o','--output_fp',help='the output filepath'),\
]

script_info['optional_options'] = [make_option('-n','--num_iterations',help='the number of iterations to perform',default=100,type='int'), make_option('-s','--sample_id_map_fp', help='Map of original sample ids to new sample ids [default: %default]', default=None)
]
script_info['version'] = __version__

comment = """# Number of entries refers to the number of rows (or cols) 
# retained in each distance matrix after filtering the distance matrices 
# to include only those samples that were in both distance matrices. 
# p-value contains the correct number of significant digits.
"""

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    sample_id_map_fp = opts.sample_id_map_fp
    if sample_id_map_fp:
        sample_id_map = dict([(k,v[0]) \
        for k,v in fields_to_dict(open(sample_id_map_fp, "U")).items()])
    else:
        sample_id_map = None

    input_dm_fps = opts.input_dms.split(',')

    output_f = open(opts.output_fp,'w')
    output_f.write(comment)
    #output_f.write('DM1\tDM2\tNumber of entries\tMantel p-value\n')

    #this is where the headin information is added, it accounts for the spacing between file names for the first two elements DM1 and DM2, but it doesn't fix the spacing between the actual number values
    output_f.write('DM1')

    num_of_first_spaces = len(input_dm_fps[0])
    first_space_string = ""
    while(num_of_first_spaces > 0):
        first_space_string = first_space_string + " "
        num_of_first_spaces = num_of_first_spaces - 1
    output_f.write(first_space_string)

    output_f.write('\tDM2')

    num_of_second_spaces = len(input_dm_fps[1])
    second_space_string = ""
    while(num_of_second_spaces > 0):
        second_space_string = second_space_string + " "
        num_of_second_spaces = num_of_second_spaces - 1
    output_f.write(second_space_string)

    num_of_entries_column_header = "Number of entries"
    output_f.write("\t" + num_of_entries_column_header)
    output_f.write('\t')

    output_f.write('Mantel p-value\n')

    num_iterations = opts.num_iterations

    for i,fp1 in enumerate(input_dm_fps):
        for fp2 in input_dm_fps[i+1:]:
            (dm1_labels, dm1), (dm2_labels, dm2) =\
             make_compatible_distance_matrices(parse_distmat(open(fp1,'U')), parse_distmat(open(fp2,'U')), lookup=sample_id_map)
            if len(dm1_labels) < 2:
                output_f.write('%s\t%s\t%d\tToo few samples\n' % (fp1,fp2,len(dm1_labels)))
                continue

            #m = Mantel(dm1, dm2, num_iterations)
            #probably better to just pass in the fp1 and fp2 here in the constructor
            m = Mantel(dm1, dm2, num_iterations)

    resultsDict = {}

    p = m.runAnalysis()

    p_str = format_p_value_for_num_iters(p,num_iterations)
    output_str = ('%s\t%s\t%d\t%s\n' % (input_dm_fps[0], input_dm_fps[1], len(dm1_labels),p_str))
    #resultsDict['Results':('%s\t%s\t%d\t%s\n' % (fp1, fp2, len(dm1_labels),p_str))]
    resultsDict['Results'] = output_str

    for line in  resultsDict['Results']:
        #easy way to output the entries
        #output_f.write(m_output[line])

        #hard way to output entries on a per item basis, modify this to alter spacing for values
        #I don't know why I didn't just use the above line because it's easier, just it doesn't handle white spacing like it should, so I wrote all the code below to handle the case for the third column so that the number is always space correctly creating a better sense of readability with the output.
        items = line.split("\t")
        tracker = 0
        for item in items:
            #This was grabbed from stack overflow in order to help with formatting...
            #http://stackoverflow.com/questions/4703390/how-to-extract-a-floating-number-from-a-string-in-python
            reFind = re.findall(r"[-+]?\d*\.\d+|\d+", item)
            if(len(reFind) > 0):
                spaces_to_add = 0
                if(tracker == 2):
                    spaces_to_add = len(num_of_entries_column_header) - len(item)
                    spaces = ""
                    if(spaces_to_add > 0):
                        while(spaces_to_add > 0):
                            spaces = spaces + " "
                            spaces_to_add = spaces_to_add - 1
                    output_f.write(item)
                    output_f.write(spaces)
                else:
                    output_f.write(item)
            else:
                output_f.write(item)

            if(tracker < len(items)-1):
                output_f.write("\t")

            tracker = tracker + 1

    output_f.close()

if __name__ == "__main__":
    main()

