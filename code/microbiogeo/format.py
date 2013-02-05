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
