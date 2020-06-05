import sys, os, unittest
import numpy as np

tau_functions_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions

class TestReadInterfaceFile(unittest.TestCase):

    def setUp(self):
        self.interface_filename = 'dummy_interface_file.dat'

        # Set reference
        self.reference_elem_connectivities = np.array([1, 2, 5, 6, 2, 3, 4, 5], dtype=int)

        # Create dummy interface file
        with open(self.interface_filename, 'w') as interface_file:
            self.WriteHeader(interface_file)
            self.WriteElementConnectivities(interface_file)


    def test_ReadInterfaceFile(self):
        # Read interface file
        with open(self.interface_filename, 'r') as interface_file:
            line = interface_file.readline()
            position_info, mesh_info, line = TauFunctions.ReadHeader(interface_file, line)
            elem_connectivities = TauFunctions.ReadElementConnectivities(interface_file, line, ElemsNr=mesh_info[1])

        # Check
        np.testing.assert_almost_equal(elem_connectivities, self.reference_elem_connectivities, decimal=16)
        self.assertIsInstance(elem_connectivities[0], int)

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


    def WriteElementConnectivities(self, interface_file):
        for i in range(len(self.reference_elem_connectivities)):
            interface_file.write(str(self.reference_elem_connectivities[i]))
            if (i + 1) % 4 == 0:
                interface_file.write('\n')
            else:
                interface_file.write(' ')


if __name__ == '__main__':
    unittest.main()