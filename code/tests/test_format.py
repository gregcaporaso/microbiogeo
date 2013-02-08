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

from microbiogeo.format import (format_method_comparison_table,
                                format_p_value_as_asterisk)

class FormatTests(TestCase):
    """Tests for the format.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.per_method_results1 = {
            'adonis': {
                'whole_body': {
                    'BODY_SITE': {
                        'full': (0.27, [0.01, 0.001]),
                        'shuffled': (0.02, [0.45, 0.476]),
                        'subsampled': ([0.24, 0.20], [[0.03, 0.005],
                                                      [0.02, 0.023]])
                    }
                }
            },

            'anosim': {
                'whole_body': {
                    'BODY_SITE': {}
                }
            }
        }

        self.per_method_results2 = {
            'adonis': self.per_method_results1['adonis'],
            'anosim': {
                'whole_body': {
                    'SEX': {}
                }
            }
        }

        self.per_method_results3 = {
            'adonis': self.per_method_results1['adonis'],
            'anosim': {
                '88_soils': {
                    'BODY_SITE': {}
                }
            }
        }

    def test_format_method_comparison_table(self):
        """Test formatting a methods summary table."""
        obs = format_method_comparison_table(self.per_method_results1)
        self.assertEqual(obs, exp_method_comparison_table)

        self.assertRaises(ValueError, format_method_comparison_table,
                          self.per_method_results2)

        self.assertRaises(ValueError, format_method_comparison_table,
                          self.per_method_results3)

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


exp_method_comparison_table = [['Method', 'whole_body\rBODY_SITE',
                                'whole_body\rBODY_SITE (shuffled)',
                                'whole_body\rBODY_SITE (subsampled)'],
                               ['adonis', '0.27; ***, ****', '0.02; x, x',
                                '0.24; **, ***\r0.20; **, **'],
                               ['anosim', ['N/A', 'N/A', 'N/A']]]


if __name__ == "__main__":
    main()
