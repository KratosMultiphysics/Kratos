from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import os, subprocess

class TestTauFunctions(KratosUnittest.TestCase):

    def test_tau_functions(self):
        # Define test file
        test_file = 'test_tau_functions_impl.py'

        # Define command for subprocess
        full_command = ['python']
        full_command.extend([test_file])

        working_directory = os.path.dirname(os.path.abspath(__file__))

        # Run test
        p = subprocess.Popen(full_command, cwd=working_directory)
        p.communicate()

if __name__ == '__main__':
    KratosUnittest.main()