import sys, os, unittest
import numpy as np

tau_functions_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions

class TestReadNodalData(unittest.TestCase):

    def setUp(self):
        # Defome dummy interface file name
        self.interface_filename = 'dummy_header_file.dat'

        # Set reference
        self.setReferenceNodalData()

        # Create dummy interface file
        with open(self.interface_filename, 'w') as interface_file:
            self.WriteNodalData(interface_file)


    def test_ReadNodalData(self):
        # Read interface file
        with open(self.interface_filename, 'r') as interface_file:
            line = interface_file.readline()
            nodal_data, line = TauFunctions.ReadNodalData(interface_file, line)

        # Check
        self.assertListEqual(nodal_data, self.reference_nodal_data)

        # Remove interface file
        os.remove(self.interface_filename)


    def setReferenceNodalData(self):
        reference_X = [0.0, 1.0, 2.0, 2.0, 1.0, 0.0]
        reference_Y = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
        reference_Z = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        reference_density = [1.2, 1.2, 1.2, 1.2, 1.2, 1.2]
        reference_pressure = [3.4, 3.4, 3.4, 3.4, 3.4, 3.4]
        self.reference_nodal_data = reference_X
        self.reference_nodal_data.extend(reference_Y)
        self.reference_nodal_data.extend(reference_Z)
        self.reference_nodal_data.extend(reference_density)
        self.reference_nodal_data.extend(reference_pressure)


    def WriteNodalData(self, interface_file):
        for i in range(len(self.reference_nodal_data)):
            interface_file.write("%.9f" % self.reference_nodal_data[i] + 'E-00')
            if (i + 1) % 5 == 0:
                interface_file.write('\n')
            else:
                interface_file.write(' ')


if __name__ == '__main__':
    unittest.main()