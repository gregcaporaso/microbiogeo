#!/usr/bin/env python

from __future__ import division

__author__ = "Damien Michael Coy"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Damien Michael Coy"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Damien Michael Coy"
__email__ = "Damien.Coy@nau.edu"
__status__ = "Development"
 
from qiime.format import format_distance_matrix
from qiime.parse import parse_distmat
from qiime.util import parse_command_line_parameters, make_option
import random

script_info = {}
script_info['brief_description'] = "Creates a distance matrix and corresponding map file"
script_info['script_description'] = """
This script creates a distance matrix and corresponding map file by taking a range \
for similar and dissimilar categories. It then creates a random number derived from the \
ranges which is stored in the distance matrix for each element. 
"""
script_info['script_usage'] = [("Create simulated data",
    "This example shows how to generate simulated data.",
    "%prog -n 8 -s 0.05,0.2 -d 0.8,1.0 -o simu_dm.txt,simu_map.txt")]
script_info['output_description'] = """
The output is the distance matrix createed by the given ranges and \
the corresponding map file.
"""
script_info['required_options'] = [
    make_option('-n','--samples', help='number of samples', type='int'),
    make_option('-s','--similar_range', help='similar lower bound, similar upper bound', type='string'),
    make_option('-d','--dissimilar_range', help='dissimilar lower bound, dissimilar upper bound', type='string'),
    make_option('-o', '--output_files', help='path to store the output distance matrix, path to store the output map', type='string')]
script_info['optional_options'] = [
    make_option('-l', '--samples_label', help='label for the samples', default="S", type='string'),
    make_option('-c', '--map_category', help='category name', default="SEX", type='string'),
    make_option('-m', '--map_labels', help='label for elements in category', default="M,F", type='string'),
    make_option('-z', '--delimiter', help='specify delimiter', default="\t", type='string')]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # Assign the input to variables
    samples = opts.samples
    similar_range = opts.similar_range.split(',')
    dissimilar_range = opts.dissimilar_range.split(',')
    output_files = opts.output_files.split(',')
    output_dm = output_files[0]
    output_map = output_files[1]
    samples_label = opts.samples_label
    map_category = opts.map_category
    map_labels = opts.map_labels.split(',')
    delimiter = opts.delimiter
    
    #setup the array of samples
    groups = [0] * samples
    index = int(samples / 2)
    while index < samples:
        groups[index] = 1
        index = index + 1
    
    #create strings that represents the distance matrix and map
    distance_matrix = ""
    distance_matrix_header = ""
    full_distance_matrix = {}
    map = ""
    map_header = "#" + samples_label + delimiter + map_category
    count = 0
    sub_count = 0
    diff = 0
    for sample in groups:
        distance_matrix_header = distance_matrix_header + delimiter + samples_label + "%s" % count
        distance_matrix = distance_matrix + "\n" + samples_label + "%s" % count
        map = map + "\n" + samples_label + "%s" % count + delimiter + map_labels[sample]
        distance_matrix_row = {} 
        while diff < count:
            distance_matrix = distance_matrix + delimiter + "%s" % full_distance_matrix[diff][count]
            diff = diff + 1
        while sub_count < samples:
            rand_num = 0
            if(sub_count == count):
                rand_num = 0
            elif(sample == groups[sub_count]):
                rand_num = round(random.uniform(float(similar_range[0]),float(similar_range[1])), 2)
            else:
                rand_num = round(random.uniform(float(dissimilar_range[0]),float(dissimilar_range[1])), 2)
            distance_matrix = distance_matrix + delimiter + "%s" % rand_num
            distance_matrix_row[sub_count] = rand_num
            sub_count = sub_count + 1
        full_distance_matrix[count] = distance_matrix_row
        count = count + 1
        sub_count = count
        diff = 0

    #write the distance matrix string to file
    try:
        file = open(output_dm, 'w')
        file.write(distance_matrix_header)
        file.write(distance_matrix)
    except:
        pass
    finally:
        file.close
    
    #write the map string to file
    try:
        file = open(output_map, 'w')
        file.write(map_header)
        file.write(map)
    except:
        pass
    finally:
        file.close

if __name__ == "__main__":
    main()

