import sys, os, unittest
import numpy as np

tau_functions_path = os.path.join(os.path.dirname(__file__), '../../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions

# Please run this file with python2. This file contains
# the implementation of the tau functions' tests

class TestTauFunctionsImpl(unittest.TestCase):

    def test_ReadInterfaceFile(self):
        self.setReference()
        self.createDummyInterfaceFile()

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

    def test_ReadElementConnectivities(self):
        number_of_elements = 2
        self.reference_elem_connectivities = np.array([1, 2, 5, 6, 2, 3, 4, 5], dtype=int)

        # Create dummy interface file
        self.interface_filename = 'dummy_interface_file.dat'
        with open(self.interface_filename, 'w') as interface_file:
            self.WriteElementConnectivities(interface_file)

        # Read interface file
        with open(self.interface_filename, 'r') as interface_file:
            line = interface_file.readline()
            elem_connectivities = TauFunctions.ReadElementConnectivities(interface_file, line, number_of_elements)

        # Check
        np.testing.assert_almost_equal(elem_connectivities, self.reference_elem_connectivities, decimal=16)
        self.assertIsInstance(elem_connectivities[0], int)

        # Remove interface file
        os.remove(self.interface_filename)

    def setReference(self):
        self.reference_position_info = ['VARIABLES', '=', '"x"', '"y"', '"z"', '"density"', '"pressure"']
        self.reference_mesh_info = [6, 2]
        reference_X = [0.0, 1.0, 2.0, 2.0, 1.0, 0.0]
        reference_Y = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
        reference_Z = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        reference_density = [1.2, 1.2, 1.2, 1.2, 1.2, 1.2]
        reference_pressure = [3.4, 3.4, 3.4, 3.4, 3.4, 3.4]
        self.reference_nodal_data = []
        self.reference_nodal_data.extend(reference_X)
        self.reference_nodal_data.extend(reference_Y)
        self.reference_nodal_data.extend(reference_Z)
        self.reference_nodal_data.extend(reference_density)
        self.reference_nodal_data.extend(reference_pressure)
        self.reference_elem_connectivities = np.array([1, 2, 5, 6, 2, 3, 4, 5], dtype=int)


    def createDummyInterfaceFile(self):
        self.interface_filename = 'dummy_interface_file.dat'
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

    def test_CalculateCellPressure(self):
        # Define nodal pressure and node ids
        nodal_pressures = [1.4, 3.6, 5.7, 7.9]
        node_ids = np.array([0, 1, 2, 3], dtype=int)

        # Compute the normal vector
        cell_pressure = TauFunctions.CalculateCellPressure(nodal_pressures, node_ids)

        # Define reference normal
        reference_pressure = 4.65

        # Check
        np.testing.assert_almost_equal(cell_pressure, reference_pressure, decimal=16)

    def test_CalculateCellArea(self):
        X, Y, Z, node_ids = self.create_dummy_cell()

        # Compute the normal vector
        cell_area = TauFunctions.CalculateCellArea(X, Y, Z, node_ids)

        # Define reference normal
        reference_area = 1.4142135623730951

        # Check
        np.testing.assert_almost_equal(cell_area, reference_area, decimal=16)

    def test_CalculateCellNormal(self):
        X, Y, Z, node_ids = self.create_dummy_cell()

        # Compute the normal vector
        cell_normal = TauFunctions.CalculateCellNormal(X, Y, Z, node_ids)

        # Define reference normal
        reference_normal = np.array([-0.7071067811865475, 0.0, 0.7071067811865475])

        # Check
        np.testing.assert_almost_equal(cell_normal, reference_normal, decimal=16)

    def create_dummy_cell(self):
        # Define cell geometry
        X = [0, 1, 1, 0]
        Y = [0, 0, 1, 1]
        Z = [0, 1, 1, 0]

        node_ids = np.array([0, 1, 2, 3], dtype=int)

        return X, Y, Z, node_ids

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

    def test_FindFilename(self):
        # Create dummy file
        path = os.getcwd() + '/'
        pattern = 'airfoilSol.pval.unsteady_i='
        step = 304
        reference_file_name = path + pattern + str(step)
        open(reference_file_name, 'w').close()

        # Retrive the file
        file_name = TauFunctions.FindFilename(path, pattern, step)

        # Check
        self.assertMultiLineEqual(file_name, reference_file_name)

        # Remove dummy file
        os.remove(file_name)

    def test_CalculateDistanceVector(self):
        # Populate two nodes with random coordinates
        X = [1.6, 19.1]
        Y = [-0.3, 3.6]
        Z = [4.36, -6.78]

        # Compute the distance vector
        distance = TauFunctions.CalculateDistanceVector(X, Y, Z, 0, 1)

        # Define reference distance
        reference_distance = np.array([17.5, 3.9, -11.14])

        # Check
        np.testing.assert_almost_equal(reference_distance, distance, decimal=16)



if __name__ == '__main__':
    unittest.main()

