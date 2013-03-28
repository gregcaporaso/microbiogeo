#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os.path import basename, join
from qiime.util import create_dir, parse_command_line_parameters, make_option

from microbiogeo.simulate import choose_gradient_subset

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = []
script_info['output_description'] = ""
script_info['required_options'] = [
    make_option('-i','--otu_table_fp', type='existing_filepath', help=''),
    make_option('-m','--map_fp', type='existing_filepath', help=''),
    make_option('-c', '--category', type='string', help=''),
    make_option('-n', '--num_total_samples', type='int', help=''),
    make_option('-o','--output_dir', type='new_dirpath', help='')
]
script_info['optional_options'] = []
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    out_dir = opts.output_dir
    create_dir(out_dir)

    subset_otu_table, subset_map_str = choose_gradient_subset(
            open(opts.otu_table_fp, 'U'), open(opts.map_fp, 'U'),
            opts.category, opts.num_total_samples)

    subset_otu_table_fp = join(out_dir, basename(opts.otu_table_fp))
    subset_otu_table_f = open(subset_otu_table_fp, 'w')
    subset_otu_table.getBiomFormatJsonString('choose_gradient_subset.py '
                                             '(microbiogeo)',
                                             subset_otu_table_f)
    subset_otu_table_f.close()

    subset_map_fp = join(out_dir, basename(opts.map_fp))
    subset_map_f = open(subset_map_fp, 'w')
    subset_map_f.write(subset_map_str)
    subset_map_f.close()


if __name__ == "__main__":
    main()
