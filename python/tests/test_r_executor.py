#!/usr/bin/env python

__author__ = "Damien Coy"
__copyright__ = "Copyright 2011, The QIIME Project" 
#remember to add yourself if you make changes
__credits__ = ["Damien Coy"] 
__license__ = "GPL"
__version__ = "1.4.0"
__maintainer__ = "Damien Coy"
__email__ = "damien.coy@nau.edu"
__status__ = "Release"


from os import remove, system, mkdir
from shutil import rmtree
from os.path import join, exists
from tempfile import NamedTemporaryFile, mkdtemp
from cogent.util.unit_test import TestCase, main
from cogent.app.util import ApplicationError
from qiime.util import get_tmp_filename

from cogent.util.misc import remove_files
from qiime.r_executor import RExecutor
from numpy import array

from test_stats import TestHelper

class RExecutorTests(TestHelper):
    """Tests of the RExecutor class"""

    def testOutput(self):
        
        # Temporary input file
        self.tmp_dm_filepath = get_tmp_filename(
            prefix='R_test_distance_matrix_',
            suffix='.txt'
            )
        seq_file = open(self.tmp_dm_filepath, 'w')

        for line in self.overview_dm_str:
            seq_file.write(line + "\n")
        seq_file.close()


        self.tmp_map_filepath = get_tmp_filename(
            prefix='R_test_map_',
            suffix='.txt'
            )
        
        seq_file = open(self.tmp_map_filepath, 'w')
        for line in self.overview_map_str:
            seq_file.write(line + "\n")
        seq_file.close()
        seq_file.close()


        self.files_to_remove = \
         [self.tmp_dm_filepath, self.tmp_map_filepath]
   
        # Prep input files in R format
        output_dir = mkdtemp()
        #self.dirs_to_remove = [output_dir]
	args = ["-d " + self.tmp_dm_filepath + " -m " + self.tmp_map_filepath + " -c DOB -o " + output_dir]

        mkdir(join(output_dir, 'rex_test'))
	rex = RExecutor()
        results = rex(args, "betadisper.r", output_dir)

#run unit tests if run from command-line
if __name__ == '__main__':
    main()
