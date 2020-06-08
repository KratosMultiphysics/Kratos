import sys, os, unittest
import numpy as np

tau_functions_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions


class TestFindInterfaceFilename(unittest.TestCase):

    def setUp(self):
        self.createTestFiles()

    def test_FindInterfaceFilename(self):
        # Retrive the file
        file_name = TauFunctions.FindInterfaceFilename(self.path, self.step)

        # Check
        self.assertMultiLineEqual(file_name, self.reference_interface_file_name)

        # Remove dummy file and directory
        TauFunctions.RemoveFilesFromPreviousSimulations()
        os.rmdir('Outputs')

    def createTestFiles(self):
        # Remove files if they exist
        TauFunctions.RemoveFilesFromPreviousSimulations()
        path = os.getcwd() + '/'
        if os.path.exists(path + "Outputs"):
            os.rmdir(path + "Outputs")
        os.mkdir('Outputs')

        # Define file names
        self.path = os.getcwd() + '/'
        pattern = 'airfoilSol.pval.unsteady_i='
        self.step = 304
        file_name = path + 'Outputs/' + pattern + str(self.step + 1)
        self.reference_interface_file_name = file_name.replace('pval', 'surface.pval') + '.dat'

        # Create dummy files
        open(file_name, 'w').close()
        open(self.reference_interface_file_name, 'w').close()


if __name__ == '__main__':
    unittest.main()

