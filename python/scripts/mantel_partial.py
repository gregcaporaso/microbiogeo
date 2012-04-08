#!/usr/bin/env python
from __future__ import division

__author__ = "Michael Dwan"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Michael Dwan"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Jai Rideout"
__email__ = "mdwan.tgen@gmail.com"
__status__ = "Development"

from os import path
from cogent.util.misc import create_dir
from qiime.parse import parse_distmat, fields_to_dict
from qiime.util import (get_options_lookup, make_compatible_distance_matrices,
                        make_option, parse_command_line_parameters)
from python.qiime.parse import DistanceMatrix
from python.qiime.stats import PartialMantel

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """
Computes a partial Mantel statistic between two distance matrices controlling for a third matrix.
"""
script_info['script_description'] = """
This script computes and outputs the p-value and the Mantel test statistic(r-value).
"""
script_info['script_usage'] = [("Compute partial Mantel statistic",
"This example computes a partial Mantel p-value on two unifrac distance matrices holding a" 
"third pH distance matrix constant and using 999 permutations in each Mantel test. Output "
"is written to the mantel_out directory.",
"%prog -i weighted_unifrac_dm.txt,unweighted_unifrac_dm.txt -c PH_dm.txt -o mantel_out "
"-n 999")]
script_info['output_description']= """
One file is created in the output directory.
"""
script_info['required_options'] = [
    make_option('-i', '--input_dms',
        help='the two input distance matrices, comma-separated'),
    make_option('-c', '--control_dm',
        help='the control matrix'),
    options_lookup['output_dir']
]
script_info['optional_options'] = [
    make_option('-n', '--num_permutations',
        help='the number of permutations to perform', default=99, type='int'),
    make_option('-s', '--sample_id_map_fp',
        help='the map of original sample ids to the new sample ids '
             '[default: %default]', default=None)
]
script_info['version'] = __version__

comment = """# p-value contains the correct number of significant digits."""

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    num_perms = opts.num_permutations

    # Try to creat the specified output dir (if not already in existence)
    try:
        create_dir(opts.output_dir)
    except:
        option_parser.error("Could not create directory specified with the -o output option.")

    input_dm_fps = opts.input_dms.split(',')
    control_dm_fp = opts.control_dm
    res_file = open(path.join(opts.output_dir, "mantel_partial_results.txt"), 'w')
    res_file.write(comment)

    # All of the matrices must be compatible for analysis to be valid.
    # credit: Greg Caporaso, from the script compare_distance_matrices.py.
    sample_id_map_fp = opts.sample_id_map_fp
    if sample_id_map_fp:
        sample_id_map = dict([(k,v[0]) \
         for k,v in fields_to_dict(open(sample_id_map_fp, "U")).items()])
    else:
        sample_id_map = None

    (dm1_labels, dm1), (dm2_labels, dm2) = make_compatible_distance_matrices(
        parse_distmat(open(input_dm_fps[0], 'U')),
        parse_distmat(open(input_dm_fps[1], 'U')), lookup=sample_id_map)

    (dm1_labels, dm1), (cdm_labels, cdm) = make_compatible_distance_matrices(
        parse_distmat(open(input_dm_fps[0], 'U')),
        parse_distmat(open(control_dm_fp, 'U')), lookup=sample_id_map)

    # Output header to result file. 
    res_file.write('\nDM1: %s\nDM2: %s\nCM: %s\npermutations: %d\n' % (input_dm_fps[0], 
                                                                     input_dm_fps[1],
                                                                     control_dm_fp,
                                                                      num_perms))

    # Construct a PartialMantel object.
    pm = PartialMantel(DistanceMatrix(dm1, dm1_labels, dm1_labels), 
                        DistanceMatrix(dm2, dm2_labels, dm2_labels), 
                        DistanceMatrix(cdm, cdm_labels, cdm_labels), num_perms)

    # Run the analysis.
    res = pm.runAnalysis()

    # Output statistic to result file.
    res_file.write('\nMantel stat(r-val)\tp-val\t')
    res_file.write('\n%f\t%f' % (res['mantel_r'], res['mantel_p']))
    res_file.close()


if __name__ == "__main__":
    main()
