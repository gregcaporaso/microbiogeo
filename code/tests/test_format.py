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

from numpy import array

from cogent.util.unit_test import TestCase, main

from microbiogeo.format import (format_method_comparison_heatmaps,
                                format_method_comparison_table)
from microbiogeo.method import Adonis, Anosim
from microbiogeo.util import StatsResults

class FormatTests(TestCase):
    """Tests for the format.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        sr_full1 = StatsResults()
        sr_full1.addResult(0.27, 0.01)
        sr_full1.addResult(0.27, 0.001)

        sr_shuff1 = StatsResults()
        sr_shuff1.addResult(0.02, 0.45)
        sr_shuff1.addResult(0.02, 0.476)

        # No empty results.
        sr_ss1 = [StatsResults(), StatsResults()]
        sr_ss1[0].addResult(0.24, 0.03)
        sr_ss1[0].addResult(0.24, 0.005)
        sr_ss1[1].addResult(0.20, 0.02)
        sr_ss1[1].addResult(0.20, 0.023)

        # Some empty results.
        sr_ss1_empty = [StatsResults(), StatsResults(), StatsResults()]
        sr_ss1_empty[0].addResult(0.24, 0.03)
        sr_ss1_empty[0].addResult(0.24, 0.005)
        sr_ss1_empty[1].addResult(0.20, 0.02)
        sr_ss1_empty[1].addResult(0.20, 0.023)

        sr_full2 = StatsResults()
        sr_full2.addResult(0.13, 0.02)
        sr_full2.addResult(0.13, 0.002)

        sr_shuff2 = StatsResults()
        sr_shuff2.addResult(0.03, 0.40)
        sr_shuff2.addResult(0.03, 0.401)

        sr_ss2 = [StatsResults(), StatsResults()]
        sr_ss2[0].addResult(0.22, 0.01)
        sr_ss2[0].addResult(0.22, 0.009)
        sr_ss2[1].addResult(0.19, 0.06)
        sr_ss2[1].addResult(0.19, 0.029)

        sr_full3 = StatsResults()
        sr_full3.addResult(0.59, 0.11)
        sr_full3.addResult(0.59, 0.101)

        sr_shuff3 = StatsResults()
        sr_shuff3.addResult(0.32, 0.65)
        sr_shuff3.addResult(0.32, 0.776)

        sr_ss3 = [StatsResults(), StatsResults()]
        sr_ss3[0].addResult(0.13, 0.23)
        sr_ss3[0].addResult(0.13, 0.105)
        sr_ss3[1].addResult(0.42, 0.92)
        sr_ss3[1].addResult(0.42, 0.723)

        sr_full4 = StatsResults()
        sr_full4.addResult(0.27, 0.01)
        sr_full4.addResult(0.27, 0.001)

        sr_shuff4 = StatsResults()
        sr_shuff4.addResult(0.02, 0.45)
        sr_shuff4.addResult(0.02, 0.476)

        sr_ss4 = [StatsResults()]
        sr_ss4[0].addResult(0.24, 0.03)
        sr_ss4[0].addResult(0.24, 0.005)

        self.per_method_results1 = {
            'adonis': {
                'whole_body': {
                    'BODY_SITE': {
                        'original': sr_full1,
                        'shuffled': sr_shuff1
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

        self.full_results1 = {
            '5_percent': {
                'unweighted_unifrac': {
                    'grouping': {
                        'adonis': {
                            'whole_body': {
                                'BODY_SITE': {
                                    'full': sr_full1,
                                    'shuffled': sr_shuff1,
                                    'subsampled': sr_ss1
                                }
                            },

                            '88_soils': {
                                'ENV_BIOME': {
                                    'full': sr_full2,
                                    'shuffled': sr_shuff2,
                                    'subsampled': sr_ss2
                                }
                            },

                            'keyboard': {}
                        },

                        'anosim': {
                            'whole_body': {
                                'BODY_SITE': {
                                    'full': sr_full1,
                                    'shuffled': sr_shuff1,
                                    'subsampled': sr_ss1
                                }
                            },

                            'keyboard': {
                                'HOST_SUBJECT_ID': {
                                    'full': sr_full3,
                                    'shuffled': sr_shuff3,
                                    'subsampled': sr_ss3
                                }
                            },

                            '88_soils': {
                                'ENV_BIOME': {}
                            }
                        }
                    }
                }
            }
        }

        # Invalid results (methods don't cover same studies).
        self.full_results2 = {
            '5_percent': {
                'unweighted_unifrac': {
                    'grouping': {
                        'adonis': {
                            'whole_body': {}
                        },

                        'anosim': {
                            '88_soils': {}
                        }
                    }
                }
            }
        }

        # Invalid results (not the same number of effect sizes to compare).
        self.full_results3 = {
            '5_percent': {
                'unweighted_unifrac': {
                    'grouping': {
                        'adonis': {
                            'whole_body': {
                                'BODY_SITE': {
                                    'full': sr_full1,
                                    'shuffled': sr_shuff1,
                                    'subsampled': sr_ss1
                                }
                            }
                        },

                        'anosim': {
                            'whole_body': {
                                'BODY_SITE': {
                                    'full': sr_full4,
                                    'shuffled': sr_shuff4,
                                    'subsampled': sr_ss4
                                }
                            }
                        }
                    }
                }
            }
        }

    def test_format_method_comparison_table(self):
        """Test formatting a methods summary table."""
        obs = format_method_comparison_table(self.per_method_results1)
        self.assertEqual(obs, exp_method_comparison_table1)

        self.assertRaises(ValueError, format_method_comparison_table,
                          self.per_method_results2)

        self.assertRaises(ValueError, format_method_comparison_table,
                          self.per_method_results3)

    def test_format_method_comparison_heatmaps(self):
        obs = format_method_comparison_heatmaps(self.full_results1,
                {'grouping': [Adonis(), Anosim()]})
        self.assertEqual(obs, exp_method_comparison_heatmaps1)

        self.assertRaises(ValueError, format_method_comparison_heatmaps,
                self.full_results2,
                {'grouping': [Adonis(), Anosim()]})

        self.assertRaises(ValueError, format_method_comparison_heatmaps,
                self.full_results3,
                {'grouping': [Adonis(), Anosim()]})


exp_method_comparison_table1 = [['Method', 'whole_body\rBODY_SITE',
                                 'whole_body\rBODY_SITE (shuffled)'],
                                ['adonis', '0.27; ***, ****', '0.02; x, x'],
                                ['anosim', 'N/A', 'N/A']]

exp_method_comparison_heatmaps1 = {
    'grouping': {
        'spearman': array([[ 1.,  1.],
                           [ 1.,  1.]]),
        'pearson': array([[ 1.,  1.],
                          [ 1.,  1.]])
    }
}


if __name__ == "__main__":
    main()
