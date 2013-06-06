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
        sr_orig1 = StatsResults()
        sr_orig1.addResult(0.27, 0.01)
        sr_orig1.addResult(0.27, 0.001)

        sr_shuff1 = StatsResults()
        sr_shuff1.addResult(0.02, 0.45)
        sr_shuff1.addResult(0.02, 0.476)

        sr_orig2 = StatsResults()
        sr_orig2.addResult(0.13, 0.02)
        sr_orig2.addResult(0.13, 0.002)

        sr_shuff2 = StatsResults()
        sr_shuff2.addResult(0.03, 0.40)
        sr_shuff2.addResult(0.03, 0.401)

        sr_orig3 = StatsResults()
        sr_orig3.addResult(0.59, 0.11)
        sr_orig3.addResult(0.59, 0.101)

        sr_shuff3 = StatsResults()
        sr_shuff3.addResult(0.32, 0.65)
        sr_shuff3.addResult(0.32, 0.776)

        sr_orig4 = StatsResults()
        sr_orig4.addResult(0.27, 0.01)
        sr_orig4.addResult(0.27, 0.001)

        sr_shuff4 = StatsResults()
        sr_shuff4.addResult(0.02, 0.45)
        sr_shuff4.addResult(0.02, 0.476)

        self.per_method_results1 = {
            'adonis': {
                'whole_body': {
                    'BODY_SITE': {
                        'original': sr_orig1,
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

        self.real_results1 = {
            '5_percent': {
                'unweighted_unifrac': {
                    'adonis': {
                        'whole_body': {
                            'BODY_SITE': {
                                'original': sr_orig1,
                                'shuffled': sr_shuff1
                            }
                        },

                        '88_soils': {
                            'ENV_BIOME': {
                                'original': sr_orig2,
                                'shuffled': sr_shuff2
                            }
                        },

                        'keyboard': {}
                    },

                    'anosim': {
                        'whole_body': {
                            'BODY_SITE': {
                                'original': sr_orig1,
                                'shuffled': sr_shuff1
                            }
                        },

                        'keyboard': {
                            'HOST_SUBJECT_ID': {
                                'original': sr_orig3,
                                'shuffled': sr_shuff3
                            }
                        },

                        '88_soils': {
                            'ENV_BIOME': {}
                        }
                    }
                }
            }
        }

        # Invalid results (methods don't cover same studies).
        self.real_results2 = {
            '5_percent': {
                'unweighted_unifrac': {
                    'adonis': {
                        'whole_body': {}
                    },

                    'anosim': {
                        '88_soils': {}
                    }
                }
            }
        }

        self.sim_results1 = {
            'adonis': {
                'whole_body': {
                    146: {
                        'BODY_SITE': {
                            1: {
                                10: {
                                    0.02: {
                                        'unweighted_unifrac': sr_orig1
                                    }
                                }
                            }
                        }
                    }
                }
            },

            'anosim': {
                'whole_body': {
                    146: {
                        'BODY_SITE': {
                            1: {
                                10: {
                                    0.02: {
                                        'unweighted_unifrac': sr_orig1
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        # Invalid results (wrong number of effect sizes).
        self.sim_results2 = {
            'adonis': {
                'whole_body': {
                    146: {
                        'BODY_SITE': {
                            1: {
                                10: {
                                    0.02: {
                                        'unweighted_unifrac': sr_orig1
                                    }
                                }
                            }
                        }
                    }
                }
            },

            'anosim': {
                'whole_body': {
                    146: {
                        'BODY_SITE': {
                            1: {
                                10: {
                                    0.02: {
                                        'unweighted_unifrac': sr_orig1,
                                        'weighted_unifrac': sr_orig2,
                                    }
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
        obs = format_method_comparison_heatmaps(self.real_results1,
                self.sim_results1, [Adonis(), Anosim()])
        self.assertEqual(obs, exp_method_comparison_heatmaps1)

        self.assertRaises(ValueError, format_method_comparison_heatmaps,
                self.real_results2, {}, [Adonis(), Anosim()])

        self.assertRaises(ValueError, format_method_comparison_heatmaps,
                self.real_results1, self.sim_results2, [Adonis(), Anosim()])


exp_method_comparison_table1 = [['Method', 'whole_body\rBODY_SITE',
                                 'whole_body\rBODY_SITE (shuffled)'],
                                ['adonis', '0.27; ***, ****', '0.02; x, x'],
                                ['anosim', 'N/A', 'N/A']]

exp_method_comparison_heatmaps1 = {
    'spearman': array([[ 1.,  1.],
                       [ 1.,  1.]]),
    'pearson': array([[ 1., 1.],
                      [ 1., 1.]])
}


if __name__ == "__main__":
    main()
