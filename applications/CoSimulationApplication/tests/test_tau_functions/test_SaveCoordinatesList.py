import sys, os, unittest
import numpy as np

tau_functions_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions

class TestSaveCoordinatesList(unittest.TestCase):

    def setUp(self):
        self.setInput()

    def test_SaveCoordinatesList(self):
        # Save coordinates list
        X, Y, Z = TauFunctions.SaveCoordinatesList(self.nodal_data, self.position_info, self.NodesNr)

        # Check
        self.assertListEqual(X, self.reference_X)
        self.assertListEqual(Y, self.reference_Y)
        self.assertListEqual(Z, self.reference_Z)

    def setInput(self):
        self.reference_X = [0.0, 1.0, 2.0, 2.0, 1.0, 0.0]
        self.reference_Y = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
        self.reference_Z = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        density = [1.2, 1.2, 1.2, 1.2, 1.2, 1.2]
        cp = [1.4, 3.4, 4.4, 7.2, 6.0, 8.1]
        self.nodal_data = []
        self.nodal_data.extend(self.reference_X)
        self.nodal_data.extend(self.reference_Y)
        self.nodal_data.extend(self.reference_Z)
        self.nodal_data.extend(density)
        self.nodal_data.extend(cp)
        self.position_info = ['VARIABLES', '=', '"x"', '"y"', '"z"', '"density"', '"cp"']
        self.NodesNr = 6
        self.velocity = 20.0

if __name__ == '__main__':
    unittest.main()