import sys, os, unittest
import numpy as np

tau_functions_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions

class TestReadInterfaceFile(unittest.TestCase):

    def setUp(self):
        self.interface_filename = 'dummy_interface_file.dat'
        self.setReference()
        self.createDummyInterfaceFile()


    def test_ReadInterfaceFile(self):
        # Read interface file
        position_info, mesh_info, nodal_data, elem_connectivities = TauFunctions.ReadInterfaceFile(self.interface_filename)

        # Check
        self.assertListEqual(position_info, self.reference_position_info)
        self.assertListEqual(mesh_info, self.reference_mesh_info)
        self.assertListEqual(nodal_data, self.reference_nodal_data)
        np.testing.assert_almost_equal(elem_connectivities, self.reference_elem_connectivities, decimal=16)
        self.assertIsInstance(elem_connectivities[0], int)

        # Remove interface file
        os.remove(self.interface_filename)


    def setReference(self):
        self.reference_position_info = ['VARIABLES', '=', '"x"', '"y"', '"z"', '"density"', '"pressure"']
        self.reference_mesh_info = [6, 2]
        self.reference_X = [0.0, 1.0, 2.0, 2.0, 1.0, 0.0]
        self.reference_Y = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
        self.reference_Z = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.reference_density = [1.2, 1.2, 1.2, 1.2, 1.2, 1.2]
        self.reference_pressure = [3.4, 3.4, 3.4, 3.4, 3.4, 3.4]
        self.reference_nodal_data = self.reference_X
        self.reference_nodal_data.extend(self.reference_Y)
        self.reference_nodal_data.extend(self.reference_Z)
        self.reference_nodal_data.extend(self.reference_density)
        self.reference_nodal_data.extend(self.reference_pressure)
        self.reference_elem_connectivities = np.array([1, 2, 5, 6, 2, 3, 4, 5], dtype=int)


    def createDummyInterfaceFile(self):
        with open(self.interface_filename, 'w') as interface_file:
            self.WriteHeader(interface_file)
            self.WriteNodalData(interface_file)
            self.WriteElementConnectivities(interface_file)


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


    def WriteNodalData(self, interface_file):
        for i in range(len(self.reference_nodal_data)):
            interface_file.write("%.9f" % self.reference_nodal_data[i] + 'E-00')
            if (i + 1) % 5 == 0:
                interface_file.write('\n')
            else:
                interface_file.write(' ')


    def WriteElementConnectivities(self, interface_file):
        for i in range(len(self.reference_elem_connectivities)):
            interface_file.write(str(self.reference_elem_connectivities[i]))
            if (i + 1) % 4 == 0:
                interface_file.write('\n')
            else:
                interface_file.write(' ')


if __name__ == '__main__':
    unittest.main()