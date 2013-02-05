#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Module to parse various supported file formats."""

# Functions to parse effect size statistics and p-values from the various
# results files.
def parse_anosim_permanova_results(results_f):
    for line in results_f:
        pass

    es, p_value = line.strip().split('\t')[1:]
    es = _parse_float(es)

    if 'Too few iters to compute p-value' in p_value:
        p_value = None
    else:
        p_value = _parse_float(p_value, 0, 1)

    return es, p_value

def parse_adonis_results(results_f):
    for line in results_f:
        if line.startswith('qiime.data$map[[opts$category]]'):
            tokens = line.strip().split()

            # The format of the file changes if the result is significant or
            # not.
            if len(tokens) == 7:
                es, p_value = line.strip().split()[-2:]
            elif len(tokens) == 8:
                es, p_value = line.strip().split()[:-1][-2:]
            else:
                raise ValueError("Encountered unparsable line: %s" % line)

            return _parse_float(es, 0, 1), _parse_float(p_value, 0, 1)

    return ValueError("Unable to parse Adonis results file.")

def parse_mrpp_results(results_f):
    a_value = None
    p_value = None

    for line in results_f:
        if line.startswith('Chance corrected within-group agreement A:'):
            tokens = line.strip().split()

            if len(tokens) == 6:
                a_value = _parse_float(tokens[-1])
            else:
                raise ValueError("Encountered unparsable line: %s" % line)
        elif line.startswith('Significance of delta:'):
            tokens = line.strip().split()

            if len(tokens) == 4:
                p_value = _parse_float(tokens[-1], 0, 1)
            else:
                raise ValueError("Encountered unparsable line: %s" % line)

    if a_value is None or p_value is None:
        raise ValueError("Unable to parse MRPP results file.")

    return a_value, p_value

def parse_dbrda_results(results_f):
    r2_value = None
    p_value = None

    for line in results_f:
        if line.startswith('Constrained'):
            tokens = line.strip().split()

            if len(tokens) == 4:
                r2_value = _parse_float(tokens[2], 0, 1)
            else:
                raise ValueError("Encountered unparsable line: %s" % line)
        elif line.startswith('Significance:'):
            tokens = line.strip().split()

            if len(tokens) == 2:
                p_value = _parse_float(tokens[1], 0, 1)
            else:
                raise ValueError("Encountered unparsable line: %s" % line)

    if r2_value is None or p_value is None:
        raise ValueError("Unable to parse db-RDA results file.")

    return r2_value, p_value

def parse_permdisp_results(results_f):
    f_value = None
    p_value = None
    at_nonparametric_section = False

    for line in results_f:
        if line.startswith('No. of permutations:'):
            at_nonparametric_section = True
        elif line.startswith('Groups') and at_nonparametric_section:
            tokens = line.strip().split()

            if len(tokens) == 7 or len(tokens) == 8:
                f_value = _parse_float(tokens[4])
                p_value = _parse_float(tokens[6], 0, 1)
            else:
                raise ValueError("Encountered unparsable line: %s" % line)

    if f_value is None or p_value is None:
        raise ValueError("Unable to parse PERMDISP results file.")

    return f_value, p_value

def parse_mantel_results(results_f):
    for line in results_f:
        pass

    if len(line.strip().split('\t')) != 7:
        raise ValueError("Encountered unparsable line: %s" % line)

    es, p_value = line.strip().split('\t')[3:5]
    es = _parse_float(es, -1, 1)
    p_value = _parse_float(p_value, 0, 1)

    return es, p_value

def parse_partial_mantel_results(results_f):
    for line in results_f:
        pass

    if len(line.strip().split('\t')) != 8:
        raise ValueError("Encountered unparsable line: %s" % line)

    es, p_value = line.strip().split('\t')[4:6]
    es = _parse_float(es, -1, 1)
    p_value = _parse_float(p_value, 0, 1)

    return es, p_value

def parse_morans_i_results(results_f):
    es_next = False
    p_value_next = False

    for line in results_f:
        if line.startswith('$observed'):
            es_next = True
        elif line.startswith('$p.value'):
            p_value_next = True
        elif es_next:
            if len(line.strip().split()) != 2:
                raise ValueError("Encountered unparsable line: %s" % line)
            es = float(line.strip().split()[1])
            if es < -1 or es > 1:
                raise ValueError("Encountered invalid I value: %.4f" % es)
            es_next = False
        elif p_value_next:
            if len(line.strip().split()) != 2:
                raise ValueError("Encountered unparsable line: %s" % line)
            p_value = float(line.strip().split()[1])
            if p_value < 0 or p_value > 1:
                raise ValueError("Encountered invalid p-value: %.4f" % p_value)
            p_value_next = False
    return es, p_value

def _parse_float(float_str, min_val=None, max_val=None):
    """Converts a float (as a string) into a float.

    Performs optional sanity checks to ensure the float is in
    [min_val, max_val].
    """
    try:
        result = float(float_str)
    except ValueError:
        raise ValueError("Could not convert float string '%s' to float."
                         % float_str)
    if (min_val is not None and result < min_val) or \
       (max_val is not None and result > max_val):
        raise ValueError("Float %.4f does not fall in valid range." % result)

    return result
