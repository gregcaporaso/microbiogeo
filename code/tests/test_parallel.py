#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Test suite for the parallel.py module."""

from os import chdir, getcwd
from os.path import exists, join
from shutil import rmtree
from tempfile import mkdtemp

from cogent.util.misc import remove_files
from cogent.util.unit_test import TestCase, main
from qiime.util import create_dir, get_qiime_temp_dir

from microbiogeo.parallel import (build_gradient_method_commands,
                                  build_gradient_method_keyboard_commands,
                                  build_grouping_method_commands)

class ParallelTests(TestCase):
    """Tests for the parallel.py module functions."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        # The prefix to use for temporary files/dirs. This prefix may be added
        # to, but all temp dirs and files created by the tests will have this
        # prefix at a minimum.
        self.prefix = 'microbiogeo_tests'

        self.start_dir = getcwd()
        self.dirs_to_remove = []
        self.files_to_remove = []

        self.tmp_dir = get_qiime_temp_dir()

        if not exists(self.tmp_dir):
            makedirs(self.tmp_dir)

            # If test creates the temp dir, also remove it.
            self.dirs_to_remove.append(self.tmp_dir)

        # Set up temporary directories to use with tests.
        self.input_dir = mkdtemp(dir=self.tmp_dir,
                                 prefix='%s_input_dir_' % self.prefix)
        self.dirs_to_remove.append(self.input_dir)

    def tearDown(self):
        """Remove temporary files/dirs created by tests."""
        # Change back to the start dir - some workflows change directory.
        chdir(self.start_dir)
        remove_files(self.files_to_remove)

        # Remove directories last, so we don't get errors trying to remove
        # files which may be in the directories.
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

    def test_build_grouping_method_commands(self):
        """Test building commands to run grouping analysis methods."""
        # Results don't exist.
        exp = [exp_grouping_method_commands1[0] % self.input_dir,
               exp_grouping_method_commands1[1] % self.input_dir]
        obs = build_grouping_method_commands(self.input_dir, '/foo/dm.txt',
                '/foo/map.txt', 'anosim', 'Treatment', [99, 999])
        self.assertEqual(obs, exp)

        # Dirs exist and are not empty.
        for permutation in 99, 999:
            tmp_results_dir = join(self.input_dir,
                                   'dm_anosim_Treatment_%d' % permutation)
            create_dir(tmp_results_dir)
            self.dirs_to_remove.append(tmp_results_dir)

            tmp_fp = join(tmp_results_dir, 'foo.txt')
            tmp_f = open(tmp_fp, 'w')
            tmp_f.write('foo')
            tmp_f.close()
            self.files_to_remove.append(tmp_fp)

        obs = build_grouping_method_commands(self.input_dir, '/foo/dm.txt',
                '/foo/map.txt', 'anosim', 'Treatment', [99, 999])
        self.assertEqual(obs, [])

        # Category dm subset that isn't the correct one.
        obs = build_grouping_method_commands(self.input_dir,
                '/foo/dm_DOB_gs3_1.txt', '/foo/map.txt', 'anosim',
                'Treatment', [99, 999])
        self.assertEqual(obs, [])

    def test_build_gradient_method_commands(self):
        """Test building commands to run gradient analysis methods."""
        # Results don't exist (mantel).
        exp = [exp_gradient_method_commands1[0] % self.input_dir,
               exp_gradient_method_commands1[1] % self.input_dir]
        obs = build_gradient_method_commands('/foo/overview', self.input_dir,
                '/foo/dm.txt', '/foo/map.txt', 'mantel', 'DOB', [999, 9999])
        self.assertEqual(obs, exp)

        # Results don't exist (morans_i).
        exp = [exp_gradient_method_commands2[0] % self.input_dir]
        obs = build_gradient_method_commands('/foo/overview', self.input_dir,
                '/foo/dm.txt', '/foo/map.txt', 'morans_i', 'DOB', [999, 9999])
        self.assertEqual(obs, exp)

        # Some dirs exist and are not empty.
        tmp_results_dir = join(self.input_dir, 'dm_mantel_DOB_9999')
        create_dir(tmp_results_dir)
        self.dirs_to_remove.append(tmp_results_dir)

        tmp_fp = join(tmp_results_dir, 'foo.txt')
        tmp_f = open(tmp_fp, 'w')
        tmp_f.write('foo')
        tmp_f.close()
        self.files_to_remove.append(tmp_fp)

        exp = [exp_gradient_method_commands1[0] % self.input_dir]
        obs = build_gradient_method_commands('/foo/overview', self.input_dir,
                '/foo/dm.txt', '/foo/map.txt', 'mantel', 'DOB', [999, 9999])
        self.assertEqual(obs, exp)

    def test_build_gradient_method_keyboard_commands(self):
        """Test building commands for running methods on keyboard study."""
        # Method doesn't apply to keyboard study.
        obs = build_gradient_method_keyboard_commands('/foo/keyboard',
                self.input_dir, '/foo/dm.txt', 'morans_i', [99, 9999])
        self.assertEqual(obs, [])

        # Results don't exist (mantel_corr).
        exp = [exp_gradient_method_keyboard_commands1[0] % self.input_dir,
               exp_gradient_method_keyboard_commands1[1] % self.input_dir]
        obs = build_gradient_method_keyboard_commands('/foo/keyboard',
                self.input_dir, '/foo/dm.txt', 'mantel_corr', [99, 9999])
        self.assertEqual(obs, exp)

        # Results don't exist (partial_mantel).
        exp = [exp_gradient_method_keyboard_commands2[0] % self.input_dir]
        obs = build_gradient_method_keyboard_commands('/foo/keyboard',
                self.input_dir, '/foo/dm.txt', 'partial_mantel', [42])
        self.assertEqual(obs, exp)


exp_grouping_method_commands1 = ['compare_categories.py --method anosim -n 99 -i /foo/dm.txt -m /foo/map.txt -c Treatment -o %s/dm_anosim_Treatment_99',
'compare_categories.py --method anosim -n 999 -i /foo/dm.txt -m /foo/map.txt -c Treatment -o %s/dm_anosim_Treatment_999']

exp_gradient_method_commands1 = ['compare_distance_matrices.py --method mantel -n 999 -i /foo/dm.txt,/foo/overview/DOB_dm.txt -o %s/dm_mantel_DOB_999',
'compare_distance_matrices.py --method mantel -n 9999 -i /foo/dm.txt,/foo/overview/DOB_dm.txt -o %s/dm_mantel_DOB_9999']

exp_gradient_method_commands2 = ['compare_categories.py --method morans_i -i /foo/dm.txt -m /foo/map.txt -c DOB -o %s/dm_morans_i_DOB']

exp_gradient_method_keyboard_commands1 = ['compare_distance_matrices.py --method mantel_corr -n 99 -i /foo/dm.txt,/foo/keyboard/euclidean_key_distances_dm.txt -o %s/dm_mantel_corr_key_distance_99',
'compare_distance_matrices.py --method mantel_corr -n 9999 -i /foo/dm.txt,/foo/keyboard/euclidean_key_distances_dm.txt -o %s/dm_mantel_corr_key_distance_9999']

exp_gradient_method_keyboard_commands2 = ['compare_distance_matrices.py --method partial_mantel -n 42 -i /foo/dm.txt,/foo/keyboard/euclidean_key_distances_dm.txt -o %s/dm_partial_mantel_key_distance_42 -c /foo/keyboard/median_unifrac_individual_distances_dm.txt']


if __name__ == "__main__":
    main()
