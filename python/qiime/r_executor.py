#!/usr/bin/env python

__author__ = "Damien Coy"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Damien Coy, Dan Knights"]
__license__ = "GPL"
__version__ = "1.4.0"
__maintainer__ = "Damien Coy"
__email__ = "damien.coy@nau.edu"
__status__ = "Release"

import subprocess
from os import remove, path, devnull
from os.path import join, abspath, dirname
from sys import stdout
from time import sleep
from tempfile import mkdtemp
from qiime.util import get_tmp_filename
from qiime.util import get_qiime_project_dir
from qiime.format import format_otu_table
from cogent.app.util import CommandLineApplication, CommandLineAppResult, \
    FilePath, ResultPath, ApplicationError
from cogent.app.parameters import Parameters
from qiime.parse import parse_otu_table, parse_mapping_file
from numpy import array, set_printoptions, nan

import sys

class RExecutor(CommandLineApplication):
    """RExecutor application controller
       Runs R with a source script (from qiime/support_files/R)
    """
    _input_handler = '_input_as_path'
    _command = "R"
    _options ={}

    _R_parameters = {
        'flags': '--slave'
        }

    # The name of the R script (located under qiime/support_files/R/)
    _R_script = ''

    _parameters = {}
    _parameters.update(_options)
    _parameters.update(_R_parameters)

    def getHelp(self):
        """Returns documentation string"""
        help_str =\
        """
        Runs the specified r script using the specified command

        Outputs:
            The results of the r script that is ran
        """
        return help_str

    def __call__(self, command_args, script_name, output_dir=None, verbose=False):
        """Run the specified r script using the commands_args

            returns a CommandLineAppResult object
        """
        input_handler = self.InputHandler
        suppress_stdout = self.SuppressStdout
        suppress_stderr = self.SuppressStderr
        if suppress_stdout:
            outfile = devnull
        else:
            outfilepath = FilePath(self.getTmpFilename(self.TmpDir))
            outfile = open(outfilepath,'w')
        if suppress_stderr:
            errfile = devnull
        else:
            errfilepath = FilePath(self.getTmpFilename(self.TmpDir))
            errfile = open(errfilepath, 'w')

	self._R_script = script_name
        rscript = self._get_R_script_path()
        base_command = self._get_base_command()
        cd_command, base_command = base_command.split(';')
        cd_command += ';'
        R_source_dir = self._get_R_script_dir()

        # Build up the command, consisting of a BaseCommand followed by
        # input and output (file) specifications

        command = self._commandline_join(
            [   cd_command, base_command,
                '--args'
            ] + command_args + [' < %s ' %(rscript)]
            )

        if self.HaltExec:
            raise AssertionError, "Halted exec with command:\n" + command

        # run command, wait for output, get exit status
        proc = subprocess.Popen(command, shell=True, stdout=outfile, stderr=errfile)

        proc.wait()
        exit_status = proc.returncode

        # Determine if error should be raised due to exit status of
        # appliciation
        if not self._accept_exit_status(exit_status):
            if exit_status == 2:
                raise ApplicationError, \
                    'R library not installed: \n' + \
                    ''.join(open(errfilepath,'r').readlines()) + '\n'
            else:
                raise ApplicationError, \
                    'Unacceptable application exit status: %s, command: %s'\
                    % (str(exit_status),command) +\
                    ' Program output: \n\n%s\n'\
                     %(''.join(open(errfilepath,'r').readlines()))

        # open the stdout and stderr if not being suppressed
        out = None
        if not suppress_stdout:
            out = open(outfilepath,"r")
        err = None
        if not suppress_stderr:
            err = open(errfilepath,"r")

	if verbose:
            msg = '\n\nCommand Executed: %s'\
                % (command) +\
                 ' \n\nR Command Output:\n%s'\
                 %(''.join(open(errfilepath,'r').readlines()))
            print(msg)

    # The methods below were taken from supervised_learning.py

    def _get_result_paths(self, output_dir):
        """Returns the filepaths for all result files"""
        files = {
            'features': 'feature_importance_scores.txt',
            'summary': 'summary.txt',
            'cv_probabilities': 'cv_probabilities.txt',
            'mislabeling': 'mislabeling.txt',
            'confusion_matrix': 'confusion_matrix.txt',
        }
        result_paths = {}
        for name, file in files.iteritems():
            result_paths[name] = ResultPath(
                Path=path.join(output_dir, file), IsWritten=True)
        return result_paths

    def _get_R_script_dir(self):
        """Returns the path to the qiime R source directory
        """

        # script_dir = path.join(qiime_dir,'qiime','support_files','R')
        # Dwan ADDED the next three lines.
        current_file_path = abspath(__file__)
        current_dir_path = dirname(current_file_path)
        script_dir = path.join(dirname(current_dir_path), '..', 'r')
        # The next two lines were the originals that Dwan replaced.
        # qiime_dir = get_qiime_project_dir()
        # script_dir = path.join(qiime_dir,'qiime','support_files','R')

        return script_dir

    def _get_R_script_path(self):
        """Returns the path to the R script to be executed
        """
        return path.join(self._get_R_script_dir(), self._R_script)

    def _commandline_join(self, tokens):
        """Formats a list of tokens as a shell command
        """
        commands = filter(None, map(str, tokens))
        return self._command_delimiter.join(commands).strip()

    def _accept_exit_status(self,exit_status):
        """ Return False to raise an error due to exit_status !=0 of application
        """
        if exit_status != 0:
            return False
        return True

    @property
    def RParameters(self):
        return self.__extract_parameters('R')

    def __extract_parameters(self, name):
        """Extracts parameters in self._<name>_parameters from self.Parameters

        Allows the program to conveniently access a subset of user-
        adjusted parameters, which are stored in the Parameters
        attribute.

        Relies on the convention of providing dicts named according to
        "_<name>_parameters" and "_<name>_synonyms".  The main
        parameters object is expected to be initialized with the
        contents of these dicts.  This method will throw an exception
        if either convention is not adhered to.
        """
        parameters = getattr(self, '_' + name + '_parameters')
        result = Parameters(parameters)
        for key in result.keys():
            result[key] = self.Parameters[key]
        return result

