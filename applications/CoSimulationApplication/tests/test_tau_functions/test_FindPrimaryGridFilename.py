import sys, os, unittest
import numpy as np

tau_functions_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions


class TestFindPrimaryGridFilename(unittest.TestCase):

    def test_FindPrimaryGridFilename(self):
        # Create dummy file
        os.mkdir('Mesh')
        path = os.getcwd() + '/'
        pattern = 'airfoil_Structured_scaliert.grid'
        reference_file_name = path + 'Mesh/' + pattern
        open(reference_file_name, 'w').close()

        # Retrive the file
        step = 200
        start_step = step
        file_name = TauFunctions.FindPrimaryGridFilename(path, step, start_step)

        # Check
        self.assertMultiLineEqual(file_name, reference_file_name)

        # Remove dummy file and directory
        os.remove(file_name)
        os.rmdir('Mesh')


if __name__ == '__main__':
    unittest.main()

