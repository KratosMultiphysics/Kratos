import sys, os, unittest
import numpy as np

tau_functions_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions

class TestSavePressure(unittest.TestCase):

    def setUp(self):
        self.setInput()

    def test_SavePressure(self):
        # Save pressure
        nodal_pressure = TauFunctions.SavePressure(self.nodal_data, self.position_info, self.NodesNr, self.velocity)

        # Define reference nodal_pressure
        reference_nodal_pressure = np.array([336.0, 816.0, 1056.0, 1728.0, 1440.0, 1944.0])

        # Check
        np.testing.assert_almost_equal(nodal_pressure, reference_nodal_pressure, decimal=16)

    def setInput(self):
        X = [0.0, 1.0, 2.0, 2.0, 1.0, 0.0]
        Y = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
        Z = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        density = [1.2, 1.2, 1.2, 1.2, 1.2, 1.2]
        cp = [1.4, 3.4, 4.4, 7.2, 6.0, 8.1]
        self.nodal_data = X
        self.nodal_data.extend(Y)
        self.nodal_data.extend(Z)
        self.nodal_data.extend(density)
        self.nodal_data.extend(cp)
        self.position_info = ['VARIABLES', '=', '"x"', '"y"', '"z"', '"density"', '"cp"']
        self.NodesNr = 6
        self.velocity = 20.0


if __name__ == '__main__':
    unittest.main()