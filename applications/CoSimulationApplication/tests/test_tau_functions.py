from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import os, subprocess, sys
import numpy as np

tau_functions_path = os.path.join(os.path.dirname(__file__), '../python_scripts/helpers')
sys.path.append(tau_functions_path)

import tau_functions as TauFunctions

class TestTauFunctions(KratosUnittest.TestCase):

    def test_RemoveFilesFromPreviousSimulations(self):
        self.createOutputsDirectory()
        os.mkdir('Mesh')

        # Create dummy empty output files
        for i in range(3):
            output_file_name = 'Outputs/test_file_' + str(i)
            open(output_file_name,'w').close()

            mesh_file_name = 'Mesh/airfoil_Structured_scaliert.grid.def' + str(i)
            open(mesh_file_name,'w').close()

        TauFunctions.RemoveFilesFromPreviousSimulations()

        # Check if files have been successfully removed
        os.rmdir('Outputs')
        os.rmdir('Mesh')


    def test_CheckIfPathExists(self):
        # Define exception
        test_filename = os.getcwd() + '/test_file.dat'
        exp_error = 'Path: "{}" not found'.format(test_filename)

        with self.assertRaisesRegex(Exception, exp_error):
            TauFunctions.CheckIfPathExists(test_filename)


    def test_WriteTautoplt(self):
        self.createTautopltTestFiles()
        self.setTautopltReferences()

        self.tautoplt_filename = TauFunctions.WriteTautoplt(self.path, self.step, self.para_path_mod, self.start_step)

        self.assertMultiLineEqual(self.tautoplt_filename, self.reference_tautoplt_filename)
        self.assertTautopltFile()

        TauFunctions.RemoveFilesFromPreviousSimulations()
        os.remove(self.reference_primary_grid_file_name)
        os.rmdir('Outputs')
        os.rmdir('Mesh')
        os.remove(self.initial_tautoplt_filename)
        os.remove(self.tautoplt_filename)


    def test_ModifyFilesIOLines(self):
        self.createTautopltTestFiles()
        self.setTautopltReferences()

        with open(self.initial_tautoplt_filename, 'r+') as tautoplt_file_reading:
            for line in tautoplt_file_reading:
                # Check if line is from IO section and test it
                line = TauFunctions.ModifyFilesIOLines(line, self.path, self.step, self.para_path_mod, self.start_step)
                if 'Primary grid filename:' in line:
                    self.assertMultiLineEqual(line, self.reference_grid_line)
                elif 'Boundary mapping filename:' in line:
                    self.assertMultiLineEqual(line, self.reference_parameter_line)
                elif 'Restart-data prefix:' in line:
                    self.assertMultiLineEqual(line, self.reference_restart_line)

            tautoplt_file_reading.close()

        TauFunctions.RemoveFilesFromPreviousSimulations()
        os.remove(self.reference_primary_grid_file_name)
        os.remove(self.initial_tautoplt_filename)
        os.rmdir('Mesh')

    def test_FindInterfaceFilename(self):
        self.createOutputsDirectory()
        self.createInterfaceFile()

        # Retrive the file
        file_name = TauFunctions.FindInterfaceFilename(self.path, self.step)

        # Check
        self.assertMultiLineEqual(file_name, self.reference_interface_file_name)

        # Remove dummy file and directory
        TauFunctions.RemoveFilesFromPreviousSimulations()
        os.rmdir('Outputs')

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
        self.assertIs(type(elem_connectivities[0]), np.int64)

        # Remove interface file
        os.remove(self.interface_filename)


    def test_SaveCoordinatesList(self):
        self.setReference()

        # Save coordinates list
        X, Y, Z = TauFunctions.SaveCoordinatesList(self.reference_nodal_data, self.reference_position_info, self.reference_mesh_info[0])

        # Check
        self.assertListEqual(X, self.reference_X)
        self.assertListEqual(Y, self.reference_Y)
        self.assertListEqual(Z, self.reference_Z)


    def test_SavePressure(self):
        self.setReference()

        # Save pressure
        nodal_pressure = TauFunctions.SavePressure(self.reference_nodal_data, self.reference_position_info, self.reference_mesh_info[0], 20.0)

        # Define reference nodal_pressure
        reference_nodal_pressure = np.array([336.0, 816.0, 1056.0, 1728.0, 1440.0, 1944.0])

        # Check
        np.testing.assert_almost_equal(nodal_pressure, reference_nodal_pressure, decimal=16)

    def test_GetCellNodeIds(self):
        # Define elem_connectivities
        elem_connectivities = np.array([24, 35, 87, 94], dtype=int)
        cell = 0

        # Compute the normal vector
        node_ids = TauFunctions.GetCellNodeIds(elem_connectivities, cell)

        # Check
        np.testing.assert_almost_equal(node_ids, elem_connectivities - 1, decimal=16)


    def test_CalculateCellForce(self):
        X, Y, Z, node_ids = self.create_dummy_cell()

        # Define nodal pressure
        nodal_pressures = [1.4, 3.6, 5.7, 7.9]

        # Compute the normal vector
        cell_force = TauFunctions.CalculateCellForce(node_ids, nodal_pressures, X, Y, Z)

        # Define reference normal
        reference_cell_force = np.array([-4.65, 0.0, 4.65])

        # Check
        np.testing.assert_almost_equal(cell_force, reference_cell_force, decimal=16)


    def test_FindPrimaryGridFilename(self):
        self.createPrimaryGridFile()

        # Retrive the file
        step = 200
        start_step = step
        file_name = TauFunctions.FindPrimaryGridFilename(self.path, step, start_step)

        # Check
        self.assertMultiLineEqual(file_name, self.reference_primary_grid_file_name)

        # Remove dummy file and directory
        os.remove(file_name)
        os.rmdir('Mesh')

    def test_FindOutputFilename(self):
        self.createOutputsDirectory()
        self.createInterfaceFile()

        # Retrive the file
        file_name = TauFunctions.FindOutputFilename(self.path, self.step)

        # Check
        self.assertMultiLineEqual(file_name, self.reference_file_name)

        # Remove dummy file and directory
        TauFunctions.RemoveFilesFromPreviousSimulations()
        os.rmdir('Outputs')


    def test_ReadHeader(self):
        # Set reference
        self.reference_position_info = ['VARIABLES', '=', '"x"', '"y"', '"z"', '"density"', '"cp"']
        self.reference_mesh_info = [6, 2]

        # Create dummy interface file
        self.interface_filename = 'dummy_header_file.dat'
        with open(self.interface_filename, 'w') as interface_file:
            self.WriteHeader(interface_file)

        # Read interface file
        with open(self.interface_filename, 'r') as interface_file:
            line = interface_file.readline()
            position_info, mesh_info, line = TauFunctions.ReadHeader(interface_file, line)

        # Check
        self.assertListEqual(position_info, self.reference_position_info)
        self.assertListEqual(mesh_info, self.reference_mesh_info)

        # Remove interface file
        os.remove(self.interface_filename)

    def test_ReadNodalData(self):
        # Set reference
        self.setReference()

        # Create dummy interface file
        self.interface_filename = 'dummy_header_file.dat'
        with open(self.interface_filename, 'w') as interface_file:
            self.WriteNodalData(interface_file)

        # Read interface file
        with open(self.interface_filename, 'r') as interface_file:
            line = interface_file.readline()
            nodal_data, line = TauFunctions.ReadNodalData(interface_file, line)

        # Check
        self.assertListEqual(nodal_data, self.reference_nodal_data)

        # Remove interface file
        os.remove(self.interface_filename)


    def test_ReadElementConnectivities(self):
        # Set reference
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
        self.assertIs(type(elem_connectivities[0]), np.int64)

        # Remove interface file
        os.remove(self.interface_filename)

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
        self.createOutputsDirectory()
        self.createInterfaceFile()

        # Retrive the file
        file_name = TauFunctions.FindFilename(self.path + 'Outputs/', self.pattern, self.step+1)

        # Check
        self.assertMultiLineEqual(file_name, self.reference_file_name)

        # Remove dummy file and directory
        TauFunctions.RemoveFilesFromPreviousSimulations()
        os.rmdir('Outputs')

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


    def createTautopltTestFiles(self):
        self.createOutputsDirectory()
        self.createInterfaceFile()
        self.createPrimaryGridFile()

        self.tautoplt_filename = self.path + 'Tautoplt.cntl'
        self.createInitialTautopltFile()
        self.step = 304
        self.para_path_mod = 'airfoil_Structured.cntl'
        self.start_step = 304


    def createInitialTautopltFile(self):
        TauFunctions.RemoveFileIfExists(self.tautoplt_filename)
        self.initial_tautoplt_filename = self.path + 'Tautoplt_initial.cntl'

        with open(self.initial_tautoplt_filename, 'w') as initial_tautoplt_file:
            initial_tautoplt_file.write('Primary grid filename:\n')
            initial_tautoplt_file.write('Boundary mapping filename:\n')
            initial_tautoplt_file.write('Restart-data prefix:\n')


    def setTautopltReferences(self):
        self.reference_tautoplt_filename = self.path + 'Tautoplt.cntl'

        primary_grid_filename = TauFunctions.FindPrimaryGridFilename(self.path, self.step, self.start_step)
        self.reference_grid_line = 'Primary grid filename:' + primary_grid_filename + ' \n'

        parameter_filename = self.path + self.para_path_mod
        self.reference_parameter_line = 'Boundary mapping filename:' + parameter_filename + ' \n'

        output_filename = TauFunctions.FindOutputFilename(self.path, self.step)
        self.reference_restart_line = 'Restart-data prefix:' + output_filename + ' \n'


    def assertTautopltFile(self):
        with open(self.tautoplt_filename, 'r+') as tautoplt_file_reading:
            # Loop over lines
            for line in tautoplt_file_reading:
                # Check if line is from IO section and modify it
                if 'Primary grid filename:' in line:
                    self.assertMultiLineEqual(line, self.reference_grid_line)
                elif 'Boundary mapping filename:' in line:
                    self.assertMultiLineEqual(line, self.reference_parameter_line)
                elif 'Restart-data prefix:' in line:
                    self.assertMultiLineEqual(line, self.reference_restart_line)

            # Close file
            tautoplt_file_reading.close()


    def setReference(self):
        self.reference_position_info = ['VARIABLES', '=', '"x"', '"y"', '"z"', '"density"', '"cp"']
        self.reference_mesh_info = [6, 2]
        self.reference_X = [0.0, 1.0, 2.0, 2.0, 1.0, 0.0]
        self.reference_Y = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
        self.reference_Z = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        reference_density = [1.2, 1.2, 1.2, 1.2, 1.2, 1.2]
        reference_cp = [1.4, 3.4, 4.4, 7.2, 6.0, 8.1]
        self.reference_nodal_data = []
        self.reference_nodal_data.extend(self.reference_X)
        self.reference_nodal_data.extend(self.reference_Y)
        self.reference_nodal_data.extend(self.reference_Z)
        self.reference_nodal_data.extend(reference_density)
        self.reference_nodal_data.extend(reference_cp)
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
        variables = 'VARIABLES = "x" "y" "z" "density" "cp"'
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


    def createOutputsDirectory(self):
        # Remove files from previous tests if they exist
        TauFunctions.RemoveFilesFromPreviousSimulations()
        path = os.getcwd() + '/'
        if os.path.exists(path + "Outputs"):
            os.rmdir(path + "Outputs")

        # Create outputs directory
        os.mkdir('Outputs')


    def createInterfaceFile(self):
        # Define file names
        self.path = os.getcwd() + '/'
        self.pattern = 'airfoilSol.pval.unsteady_i='
        self.step = 304
        self.reference_file_name = self.path + 'Outputs/' + self.pattern + str(self.step + 1)
        self.reference_interface_file_name = self.reference_file_name.replace('pval', 'surface.pval') + '.dat'

        # Create dummy files
        open(self.reference_file_name, 'w').close()
        open(self.reference_interface_file_name, 'w').close()


    def createPrimaryGridFile(self):
        # Create dummy file
        os.mkdir('Mesh')
        self.path = os.getcwd() + '/'
        pattern = 'airfoil_Structured_scaliert.grid'
        self.reference_primary_grid_file_name = self.path + 'Mesh/' + pattern
        open(self.reference_primary_grid_file_name, 'w').close()


    def create_dummy_cell(self):
        # Define cell geometry
        X = [0, 1, 1, 0]
        Y = [0, 0, 1, 1]
        Z = [0, 1, 1, 0]

        node_ids = np.array([0, 1, 2, 3], dtype=int)

        return X, Y, Z, node_ids


if __name__ == '__main__':
    KratosUnittest.main()