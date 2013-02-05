#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Test suite for the format.py module."""

from cogent.util.unit_test import TestCase, main

from microbiogeo.format import format_p_value_as_asterisk

class FormatTests(TestCase):
    """Tests for the format.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        pass

    def test_format_p_value_as_asterisk(self):
        """Test formatting a p-value to indicate statistical significance."""
        obs = format_p_value_as_asterisk(1.0)
        self.assertEqual(obs, 'x')

        obs = format_p_value_as_asterisk(0.09)
        self.assertEqual(obs, '*')

        obs = format_p_value_as_asterisk(0.045)
        self.assertEqual(obs, '**')

        obs = format_p_value_as_asterisk(0.01)
        self.assertEqual(obs, '***')

        obs = format_p_value_as_asterisk(0.0005)
        self.assertEqual(obs, '****')

    def test_format_p_value_as_asterisk_invalid_input(self):
        """Test supplying an invalid p-value results in error being thrown."""
        self.assertRaises(TypeError, format_p_value_as_asterisk, 1)
        self.assertRaises(TypeError, format_p_value_as_asterisk, "0.05")
        self.assertRaises(TypeError, format_p_value_as_asterisk, [0.05])

        self.assertRaises(ValueError, format_p_value_as_asterisk, 1.1)
        self.assertRaises(ValueError, format_p_value_as_asterisk, -0.042)


if __name__ == "__main__":
    main()
