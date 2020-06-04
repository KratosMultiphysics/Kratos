import sys, os, unittest
import numpy as np

tau_functions_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions


class TestFindInitialMeshFilename(unittest.TestCase):

    def test_FindInitialMeshFilename(self):
        # Create dummy file
        path = os.getcwd() + '/'
        pattern = 'airfoil_Structured_scaliert.grid'
        reference_file_name = path + pattern
        open(reference_file_name, 'w').close()

        # Retrive the file
        file_name = TauFunctions.FindInitialMeshFilename(path, pattern)

        # Check
        self.assertMultiLineEqual(file_name, reference_file_name)

        # Remove dummy file
        os.remove(file_name)


if __name__ == '__main__':
    unittest.main()

