from qiime.util import make_option
from parse import parse_distmat
from qiime.format import format_p_value_for_num_iters
from qiime.util import (parse_command_line_parameters, 
                        get_options_lookup,
                        make_compatible_distance_matrices)

from numpy import array, asarray, ravel, sqrt
from numpy.random import permutation

from mantel import Mantel

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
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
    output_f.write('DM1\tDM2\tNumber of entries\tMantel p-value\n')

    num_iterations = opts.num_iterations

    m = Mantel(sample_id_map, input_dm_fps, num_iterations)

    for line in m.runAnalysis():
        output_f.write(line)
    output_f.close()

if __name__ == "__main__":
    main()

