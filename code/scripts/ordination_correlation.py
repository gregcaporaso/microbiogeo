#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os.path import join
from qiime.util import create_dir, parse_command_line_parameters, make_option

from microbiogeo.ordination_correlation import (compute_ordination_correlation,
                                                CORRELATION_TYPES)

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = []
script_info['output_description'] = ""
script_info['required_options'] = [
    make_option('-i','--coord_fp', type='existing_filepath', help=''),
    make_option('-m','--map_fp', type='existing_filepath', help=''),
    make_option('-c', '--category', type='string', help='must be numeric'),
    make_option('-o','--output_dir', type='new_dirpath', help='')
]
script_info['optional_options'] = [
    make_option('-a', '--axis', type='int', help='the axis number. Must be '
        'between 1 and the number of axes in the input coordinates file '
        '(inclusive) [default: %default]', default=1),
    make_option('-t', '--type', type='choice',
        choices=CORRELATION_TYPES, help='the type of correlation coefficient '
        'to compute. Valid choices: ' + ' or '.join(CORRELATION_TYPES) +
        ' [default: %default]', default='pearson'),
    make_option('-n','--num_permutations', type='int',
        help='the number of permutations to perform when calculating the '
        'nonparametric p-value. Must be an integer greater than or equal to '
        'zero. If zero, the nonparametric test of significance will not be '
        'performed and the nonparametric p-value will be reported as "N/A" '
        '[default: %default]', default=999)
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    output_dir = opts.output_dir
    create_dir(output_dir)

    map_f = open(opts.map_fp, 'U')
    coord_f = open(opts.coord_fp, 'U')

    results = compute_ordination_correlation(map_f, coord_f, opts.category,
            opts.axis, opts.type, opts.num_permutations)

    header = ['Correlation coefficient', 'Parametric p-value',
              'Nonparametric p-value']
    with open(join(output_dir, 'ord_corr_results.txt'), 'w') as output_f:
        output_f.write('\t'.join(header) + '\n')
        output_f.write('%.4f\t%.4f\t%s\n' % results)


if __name__ == "__main__":
    main()
