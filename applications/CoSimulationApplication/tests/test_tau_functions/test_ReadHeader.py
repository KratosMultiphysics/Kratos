import sys, os, unittest
import numpy as np

tau_functions_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions

class TestReadHeader(unittest.TestCase):

    def setUp(self):
        # Defome dummy interface file name
        self.interface_filename = 'dummy_header_file.dat'

        # Set reference
        self.reference_position_info = ['VARIABLES', '=', '"x"', '"y"', '"z"', '"density"', '"pressure"']
        self.reference_mesh_info = [6, 2]

        # Create dummy interface file
        with open(self.interface_filename, 'w') as interface_file:
            self.WriteHeader(interface_file)


    def test_ReadHeader(self):
        # Read interface file
        with open(self.interface_filename, 'r') as interface_file:
            line = interface_file.readline()
            position_info, mesh_info, line = TauFunctions.ReadHeader(interface_file, line)

        # Check
        self.assertListEqual(position_info, self.reference_position_info)
        self.assertListEqual(mesh_info, self.reference_mesh_info)

        # Remove interface file
        os.remove(self.interface_filename)


    def WriteHeader(self, interface_file):
        title = 'TITLE     = '
        title += '"Grid: /path/to/Mesh/airfoil_Structured_scaliert.grid, '
        title += 'Pointdata: path/to/Outputs/airfoilSol.pval.unsteady_i=201_t=1.00500e+00"'
        variables = 'VARIABLES = "x" "y" "z" "density" "pressure"'
        zone = 'ZONE T="MEMBRANE"'
        mesh_info = 'N=6, E=2, F=FEBLOCK ET=Quadrilateral'
        interface_file_lines = [title, variables, zone, mesh_info]
        for i in range(len(interface_file_lines)):
            interface_file.write(interface_file_lines[i] + '\n')
        interface_file.write('\n')


if __name__ == '__main__':
    unittest.main()