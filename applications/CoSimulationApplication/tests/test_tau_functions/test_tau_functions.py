from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import os, subprocess

class TestTauFunctions(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def test_RemoveFilesFromPreviousSimulations(self):
        # Define test file
        test_file = 'test_RemoveFilesFromPreviousSimulations.py'

        # Run function
        self.execute_test(test_file)

    def test_ReadInterfaceFile(self):
        # Define test file
        test_file = 'test_ReadInterfaceFile.py'

        # Run test
        self.execute_test(test_file)

    def test_SavePressure(self):
        # Define test file
        test_file = 'test_SavePressure.py'

        # Run test
        self.execute_test(test_file)

    def test_GetCellNodeIds(self):
        # Define test file
        test_file = 'test_GetCellNodeIds.py'

        # Run test
        self.execute_test(test_file)

    def test_CalculateCellForce(self):
        # Define test file
        test_file = 'test_CalculateCellForce.py'

        # Run test
        self.execute_test(test_file)

    def test_FindPrimaryGridFilename(self):
        # Define test file
        test_file = 'test_FindPrimaryGridFilename.py'

        # Run test
        self.execute_test(test_file)

    def test_FindOutputFilename(self):
        # Define test file
        test_file = 'test_FindOutputFilename.py'

        # Run test
        self.execute_test(test_file)

    def test_ReadHeader(self):
        # Define test file
        test_file = 'test_ReadHeader.py'

        # Run test
        self.execute_test(test_file)

    def test_ReadNodalData(self):
        # Define test file
        test_file = 'test_ReadNodalData.py'

        # Run test
        self.execute_test(test_file)

    def test_ReadElementConnectivities(self):
        # Define test file
        test_file = 'test_ReadElementConnectivities.py'

        # Run test
        self.execute_test(test_file)

    def test_CalculateCellPressure(self):
        # Define test file
        test_file = 'test_CalculateCellPressure.py'

        # Run test
        self.execute_test(test_file)

    def test_CalculateCellArea(self):
        # Define test file
        test_file = 'test_CalculateCellArea.py'

        # Run test
        self.execute_test(test_file)

    def test_CalculateCellNormal(self):
        # Define test file
        test_file = 'test_CalculateCellNormal.py'

        # Run test
        self.execute_test(test_file)

    def test_FindInitialMeshFilename(self):
        # Define test file
        test_file = 'test_FindInitialMeshFilename.py'

        # Run test
        self.execute_test(test_file)

    def test_FindFilename(self):
        # Define test file
        test_file = 'test_FindFilename.py'

        # Run test
        self.execute_test(test_file)

    def test_CalculateDistanceVector(self):
        # Define test file
        test_file = 'test_CalculateDistanceVector.py'

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
