from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import os, subprocess
# import execnet

# import imp
# foo = imp.load_source('module.name', '/path/to/file.py')
# import KratosMultiphysics.CoSimulationApplication.helpers.tau_functions as tau_functions


class TestTauFunctions(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def test_RemoveFilesFromPreviousSimulations(self):
        # Create dummy outputs and mesh directory
        os.mkdir('Outputs')
        os.mkdir('Mesh')

        # Create dummy empty output files
        for i in range(3):
            output_file_name = 'Outputs/test_file_' + str(i)
            open(output_file_name,'w').close()

            mesh_file_name = 'Mesh/airfoil_Structured_scaliert.grid.def' + str(i)
            open(mesh_file_name,'w').close()

        # Define test file
        test_file = 'test_tau_functions/test_RemoveFilesFromPreviousSimulations.py'

        # Run function
        self.execute_test(test_file)

        # Check if files have been successfully removed
        os.rmdir('Outputs')
        os.rmdir('Mesh')

    def test_CalculateDistanceVector(self):
        # Define test file
        test_file = 'test_tau_functions/test_CalculateDistanceVector.py'

        # Run test
        self.execute_test(test_file)

    def execute_test(self, test_file):
        full_command = ['python']
        full_command.extend([test_file])

        working_directory = os.path.dirname(os.path.abspath(__file__))

        p = subprocess.Popen(full_command, cwd=working_directory)
        p.communicate()

if __name__ == '__main__':
    KratosUnittest.main()
