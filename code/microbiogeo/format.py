#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Module to format data for presentation."""

# Not unit-tested.
def create_results_summary_tables(results, out_dir, filename_prefix):
    for depth_desc, depth_res in results.items():
        for metric, metric_res in depth_res.items():
            table_rows = format_method_comparison_table(metric_res)

            with open(join(out_dir, '%s_%s_%s.txt' % (filename_prefix,
                                                      depth_desc,
                                                      metric)), 'wb') as out_f:
                # We use \r so that we can force linebreaks within cells when
                # imported into Excel. Not sure if this will work with other
                # spreadsheet programs such as Open Office.
                tsv_writer = writer(out_f, delimiter='\t', lineterminator='\r')
                tsv_writer.writerows(table_rows)

def format_method_comparison_table(per_method_results):
    constructed_header = False
    header = ['Method']
    rows = []

    for method, method_res in sorted(per_method_results.items()):
        row = [method]

        for study, study_res in sorted(method_res.items()):
            for category, category_res in sorted(study_res.items()):
                header_column_title = '%s\r%s' % (study, category)
                header.append(header_column_title)
                header.append(header_column_title + ' (shuffled)')
                header.append(header_column_title + ' (subsampled)')

                if len(category_res) > 0:
                    # Format full results.
                    row.append('%.2f; %s' % (category_res['full'][0],
                               ', '.join(map(format_p_value_as_asterisk,
                                             category_res['full'][1]))))

                    # Format shuffled results.
                    row.append('%.2f; %s' % (category_res['shuffled'][0],
                               ', '.join(map(format_p_value_as_asterisk,
                                             category_res['shuffled'][1]))))

                    # Format subsampled results.
                    cell = ['%.2f; %s' % (es, ', '.join(map(
                            format_p_value_as_asterisk, p_vals)))
                            for es, p_vals in zip(*category_res['subsampled'])]
                    row.append('\r'.join(cell))
                else:
                    row.append(['N/A'] * 3)

        if not constructed_header:
            rows.append(header[:])
            constructed_header = True

        rows.append(row)

    return rows

def format_p_value_as_asterisk(p_value):
    if not isinstance(p_value, float):
        raise TypeError("p-value must be a float.")
    if p_value < 0 or p_value > 1:
        raise ValueError("p-value must be a float between 0 and 1, inclusive.")

    result = 'x'

    if p_value <= 0.1:
        result = '*'
    if p_value <= 0.05:
        result += '*'
    if p_value <= 0.01:
        result += '*'
    if p_value <= 0.001:
        result += '*'

    return result
