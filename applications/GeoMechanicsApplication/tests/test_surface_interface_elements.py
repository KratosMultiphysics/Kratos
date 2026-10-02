import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsSurfaceInterfaceElementTests(KratosUnittest.TestCase):
    """
    This class contains one test case for a 3+3, 4+4, 6+6 and 8+8 surface interface elements which are loaded
    in compression as well as in shear. This is done both by applying a load (Neumann)
    and prescribing a displacement (Dirichlet) on one side of each interface.
    """
    def setUp(self):
        self.model = Kratos.Model()
        self.normal_stiffness = 3.0e7
        self.shear_stiffness = 1.5e7

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def assertVectorsAlmostEqual(self, actual_vectors, expected_vectors):
        for actual_vector, expected_vector in zip(actual_vectors, expected_vectors):
            self.assertVectorAlmostEqual(actual_vector, expected_vector)

    def test_neumann_single_stage_3_plus_3(self):
        simulation = self.run_simulation('Neumann_single_stage_3_plus_3')
        self.assert_outputs_for_surface_interface_element(simulation, 4, 3, 6)

    def test_neumann_single_stage_3_plus_3_umat(self):
        simulation = self.run_simulation('Neumann_single_stage_3_plus_3_umat')
        self.assert_outputs_for_surface_interface_element(simulation, 4, 3, 6)

    def test_neumann_single_stage_4_plus_4(self):
        simulation = self.run_simulation('Neumann_single_stage_4_plus_4')
        self.assert_outputs_for_surface_interface_element(simulation, 4, 4, 3)

    @KratosUnittest.skip("test skipped as results are unexpected")
    def test_neumann_single_stage_6_plus_6(self):
        simulation = self.run_simulation('Neumann_single_stage_6_plus_6')
        self.assert_outputs_for_surface_interface_element(simulation, 9, 3, 6)

    @KratosUnittest.skip("test skipped as results are unexpected")
    def test_neumann_single_stage_8_plus_8(self):
        simulation = self.run_simulation('Neumann_single_stage_8_plus_8')
        self.assert_outputs_for_surface_interface_element(simulation, 8, 4, 3)

    def test_dirichlet_single_stage_3_plus_3(self):
        simulation = self.run_simulation('Dirichlet_single_stage_3_plus_3')
        self.assert_outputs_for_surface_interface_element(simulation, 4, 3, 6)

    def test_dirichlet_single_stage_4_plus_4(self):
        simulation = self.run_simulation('Dirichlet_single_stage_4_plus_4')
        self.assert_outputs_for_surface_interface_element(simulation, 4, 4, 3)
    
    def test_dirichlet_single_stage_6_plus_6(self):
        simulation = self.run_simulation('Dirichlet_single_stage_6_plus_6')
        self.assert_outputs_for_surface_interface_element(simulation, 9, 3, 6)
    
    def test_dirichlet_single_stage_8_plus_8(self):
        simulation = self.run_simulation('Dirichlet_single_stage_8_plus_8')
        self.assert_outputs_for_surface_interface_element(simulation, 8, 4, 3)

    def test_neumann_multi_stage_3_plus_3(self):
        self.run_simulation_multistages('Neumann_multi_stage_3_plus_3', 4, 3, [-1334.0, -3000.0], [-666.0, -1332.0])

    def test_neumann_multi_stage_3_plus_3_umat(self):
        self.run_simulation_multistages('Neumann_multi_stage_3_plus_3_umat', 4, 3, [-1334.0, -3000.0], [-666.0, -1332.0])

    def test_neumann_multi_stage_4_plus_4(self):
        self.run_simulation_multistages('Neumann_multi_stage_4_plus_4', 4, 4, [-1334.0, -3000.0], [-666.0, -1332.0])

    @KratosUnittest.skip("test skipped as results are unexpected")
    def test_neumann_multi_stage_6_plus_6(self):
        self.run_simulation_multistages('Neumann_multi_stage_6_plus_6', 9, 3, [-1334.0, -3000.0], [-666.0, -1332.0])

    @KratosUnittest.skip("test skipped as results are unexpected")
    def test_neumann_multi_stage_8_plus_8(self):
        self.run_simulation_multistages('Neumann_multi_stage_8_plus_8', 8, 4, [-1334.0, -3000.0], [-666.0, -1332.0])

    def test_dirichlet_multi_stage_3_plus_3(self):
        self.run_simulation_multistages('Dirichlet_multi_stage_3_plus_3', 4, 3, [-1334.0, -3000.0], [-666.0, -1332.0])

    def test_dirichlet_multi_stage_4_plus_4(self):
        self.run_simulation_multistages('Dirichlet_multi_stage_4_plus_4', 4, 4, [-1334.0, -3000.0], [-666.0, -1332.0])

    def test_dirichlet_multi_stage_6_plus_6(self):
        self.run_simulation_multistages('Dirichlet_multi_stage_6_plus_6', 9, 3, [-1334.0, -3000.0], [-666.0, -1332.0])

    def test_dirichlet_multi_stage_8_plus_8(self):
        self.run_simulation_multistages('Dirichlet_multi_stage_8_plus_8', 8, 4, [-1334.0, -3000.0], [-666.0, -1332.0])


    @staticmethod
    def run_simulation(test_name):
        file_path = test_helper.get_file_path(os.path.join('surface_interface_elements', test_name))
        return test_helper.run_kratos(file_path)
        
    def run_simulation_multistages(self, test_name, number_of_nodes, number_of_integration_points, shear_tractions, normal_tractions):
        file_path = test_helper.get_file_path(os.path.join('surface_interface_elements', test_name))
        initial_cwd = os.getcwd()
        os.chdir(file_path)
        project_parameters_file_names = ['ProjectParameters_stage1.json', 'ProjectParameters_stage2.json']
        for file_name, normal_traction, shear_traction in zip(project_parameters_file_names, normal_tractions, shear_tractions):
            simulation = test_helper.make_geomechanics_analysis(self.model, os.path.join(file_path, file_name))
            simulation.Run()
            self.assert_outputs_for_surface_interface_element_multistages(simulation, number_of_nodes, number_of_integration_points, shear_traction, normal_traction)
        os.chdir(initial_cwd)


    def assert_outputs_for_surface_interface_element(self, simulation, number_of_nodes, number_of_integration_points, number_of_elements):
        displacements = test_helper.get_displacement(simulation)
        tractions = test_helper.get_on_integration_points(simulation, KratosGeo.GEO_EFFECTIVE_TRACTION_VECTOR)
        relative_displacements = test_helper.get_on_integration_points(simulation, KratosGeo.GEO_RELATIVE_DISPLACEMENT_VECTOR)

        shear_traction = -667.0
        expected_shear_displacement = shear_traction / self.shear_stiffness
        normal_traction = -333.0
        expected_normal_displacement = normal_traction / self.normal_stiffness

        # Check the element at XY
        # 5*number_of_nodes and 6*number_of_nodes refer to start and end of nodes for a specific element (here in XY plane)
        # Arranged in mdpa file with such specific indices
        for index in range(5*number_of_nodes, 6*number_of_nodes):
            self.assertAlmostEqual(displacements[index][0], expected_shear_displacement)
            self.assertAlmostEqual(displacements[index][1], expected_shear_displacement)
            self.assertAlmostEqual(displacements[index][2], expected_normal_displacement)

        node_index = number_of_elements * 2 // 3
        tractions_horizontal_element = tractions[node_index]
        for index in range(number_of_integration_points):
            self.assertAlmostEqual(tractions_horizontal_element[index][0], normal_traction)
            self.assertAlmostEqual(tractions_horizontal_element[index][1], shear_traction)
            self.assertAlmostEqual(tractions_horizontal_element[index][2], shear_traction)
            
        relative_displacements_horizontal_element = relative_displacements[node_index]
        for index in range(number_of_integration_points):
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][0], expected_normal_displacement)
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][1], expected_shear_displacement)
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][2], expected_shear_displacement)
  
        # Check the element at YZ
        for index in range(3*number_of_nodes, 4*number_of_nodes):
            self.assertAlmostEqual(displacements[index][0], expected_normal_displacement)
            self.assertAlmostEqual(displacements[index][1], expected_shear_displacement)
            self.assertAlmostEqual(displacements[index][2], expected_shear_displacement)
        
        node_index = number_of_elements // 3
        tractions_horizontal_element = tractions[node_index]
        for index in range(number_of_integration_points):
            self.assertAlmostEqual(tractions_horizontal_element[index][0], normal_traction)
            self.assertAlmostEqual(tractions_horizontal_element[index][1], shear_traction)
            self.assertAlmostEqual(tractions_horizontal_element[index][2], shear_traction)

        relative_displacements_horizontal_element = relative_displacements[node_index]
        for index in range(number_of_integration_points):
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][0], expected_normal_displacement)
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][1], expected_shear_displacement)
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][2], expected_shear_displacement)

        # Check the element at ZX
        for index in range(number_of_nodes, 2*number_of_nodes):
            self.assertAlmostEqual(displacements[index][0], expected_shear_displacement)
            self.assertAlmostEqual(displacements[index][1], expected_normal_displacement)
            self.assertAlmostEqual(displacements[index][2], expected_shear_displacement)

        node_index = 0
        tractions_horizontal_element = tractions[node_index]
        for index in range(number_of_integration_points):
            self.assertAlmostEqual(tractions_horizontal_element[index][0], normal_traction)
            self.assertAlmostEqual(tractions_horizontal_element[index][1], shear_traction)
            self.assertAlmostEqual(tractions_horizontal_element[index][2], shear_traction)

        relative_displacements_horizontal_element = relative_displacements[node_index]
        for index in range(number_of_integration_points):
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][0], expected_normal_displacement)
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][1], expected_shear_displacement)
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][2], expected_shear_displacement)

  
    def assert_outputs_for_surface_interface_element_multistages(self, simulation, number_of_nodes, number_of_integration_points, shear_traction, normal_traction):
        displacements = test_helper.get_displacement(simulation)

        expected_shear_displacement = shear_traction / self.shear_stiffness
        expected_normal_displacement = normal_traction / self.normal_stiffness

        # Check the element at ZX
        for index in range(number_of_nodes, 2*number_of_nodes):
            self.assertAlmostEqual(displacements[index][0], expected_shear_displacement)
            self.assertAlmostEqual(displacements[index][1], expected_normal_displacement)
            self.assertAlmostEqual(displacements[index][2], expected_shear_displacement)

        tractions = test_helper.get_on_integration_points(simulation, KratosGeo.GEO_EFFECTIVE_TRACTION_VECTOR)
        node_index = 0
        tractions_horizontal_element = tractions[node_index]
        for index in range(number_of_integration_points):
            self.assertAlmostEqual(tractions_horizontal_element[index][0], normal_traction)
            self.assertAlmostEqual(tractions_horizontal_element[index][1], shear_traction)
            self.assertAlmostEqual(tractions_horizontal_element[index][2], shear_traction)

        relative_displacements = test_helper.get_on_integration_points(simulation, KratosGeo.GEO_RELATIVE_DISPLACEMENT_VECTOR)
        relative_displacements_horizontal_element = relative_displacements[node_index]
        for index in range(number_of_integration_points):
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][0], expected_normal_displacement)
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][1], expected_shear_displacement)
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][2], expected_shear_displacement)
  
  
if __name__ == '__main__':
    KratosUnittest.main()
