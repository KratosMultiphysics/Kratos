import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsSurfaceInterfaceElementTests(KratosUnittest.TestCase):
    """
    This class contains one test case for a 3+3 line interface element which is loaded
    in compression as well as in shear. This is done both by applying a load (Neumann)
    and prescribing a displacement (Dirichlet) on one side of the interface.
    """
    def setUp(self):
        self.model = Kratos.Model()
        self.normal_stiffness = 3.0e7
        self.shear_stiffness = 1.5e7

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass


    #def check_displacement_single_stage(self, test_name, number_of_nodes):
    #    shear_traction = 667.0
    #    expected_horizontal_displacement = -shear_traction / self.shear_stiffness
    #    normal_traction = 333.0
    #    expected_vertical_displacement = -normal_traction / self.normal_stiffness
    #    displacement = test_helper.get_displacement(simulation)
    #    for base in [10, 30, 50]:
    #        for index in range(base, base + number_of_nodes):
    #            self.assertAlmostEqual(displacements[index][0], expected_horizontal_displacement)
    #            self.assertAlmostEqual(displacements[index][1], expected_vertical_displacement)
    #
    #
    #def check_tractions_single_stage(self, test_name, number_of_integration_points):
    #    normal_displacement = 1.11e-5
    #    tangential_displacement_1 = -4.4466666666666666e-5
    #    tangential_displacement_2 = -4.4466666666666666e-5
    #    expected_normal_traction = normal_displacement * self.normal_stiffnes
    #    expected_shear_traction_1 = tangential_displacement_1 * self.shear_stiffness
    #    expected_shear_traction_2 = tangential_displacement_2 * self.shear_stiffness
    #    file_path = test_helper.get_file_path(os.path.join('surface_interface_elements', test_name))
    #    simulation = test_helper.run_kratos(file_path)
    #    tractions = test_helper.get_on_integration_points(simulation, Kratos.CAUCHY_STRESS_VECTOR)
    #    for index in range(number_of_integration_points):
    #        self.assertAlmostEqual(stresses[index][0], expected_normal_traction)
    #        self.assertAlmostEqual(stresses[index][1], expected_shear_traction_1)
    #        self.assertAlmostEqual(stresses[index][2], expected_shear_traction_2)
    #
    #
    #def check_displacement_multi_stage(self, test_name, number_of_nodes):
    #    project_parameters_file_names = ['ProjectParameters_stage1.json', 'ProjectParameters_stage2.json']
    #    for file_name in (project_parameters_file_names):
    #        stage = test_helper.make_geomechanics_analysis(self.model, os.path.join('surface_interface_elements', test_name))
    #        stage.Run()
    #    traction_vectors = [[666.0, 1334.0, 1334.0], [1332.0, 3000.0, 3000.0]]
    #    expected_displacement_vectors = [[traction_vectors[0] / self.normal_stiffness,
    #                                      traction_vectors[1] / self.shear_stiffness]] * number_of_nodes
    #    self.assertVectorsAlmostEqual(test_helper.get_displacement(stage), expected_displacement_vectors)
    #
    #
    #def check_tractions_multi_stage(self, test_name, number_of_integration_points):
    #    project_parameters_file_names = ['ProjectParameters_stage1.json', 'ProjectParameters_stage2.json']
    #    for file_name in (project_parameters_file_names):
    #        stage = test_helper.make_geomechanics_analysis(self.model, os.path.join('surface_interface_elements', test_name))
    #        stage.Run()
    #    displacement_vectors = [[-8.8933333333333332e-5, -2.22e-5, -2.22e-5], [-2.0e-4, -4.44e-5, -4.44e-5]]
    #    expected_traction_vectors = [[self.normal_stiffness * displacement_vector[0],
    #                                  self.shear_stiffness * displacement_vector[1]]] * number_of_nodes
    #    self.assertVectorsAlmostEqual(test_helper.get_on_integration_points(stage, Kratos.CAUCHY_STRESS_VECTOR)[0],
    #                                  expected_traction_vectors)


    def assertVectorsAlmostEqual(self, actual_vectors, expected_vectors):
        for actual_vector, expected_vector in zip(actual_vectors, expected_vectors):
            self.assertVectorAlmostEqual(actual_vector, expected_vector)


    #def test_neumann_single_stage_3_plus_3(self):
    #    simulation = self.run_simulation('neumann_single_stage_3_plus_3')
    #    self.assert_outputs_for_surface_interface_element(simulation, 4, 3, 6)
    #
    #def test_neumann_single_stage_4_plus_4(self):
    #    self.check_displacement_single_stage('neumann_single_stage_4_plus_4', 4)
    #
    #def test_neumann_single_stage_6_plus_6(self):
    #    self.check_displacement_single_stage('neumann_single_stage_6_plus_6', 9)
    #
    #def test_neumann_single_stage_8_plus_8(self):
    #    self.check_displacement_single_stage('neumann_single_stage_8_plus_8', 8)

    def test_dirichlet_single_stage_3_plus_3(self):
        simulation = self.run_simulation('dirichlet_single_stage_3_plus_3')
        self.assert_outputs_for_surface_interface_element(simulation, 4, 3, 6)

    def test_dirichlet_single_stage_4_plus_4(self):
        simulation = self.run_simulation('dirichlet_single_stage_4_plus_4')
        self.assert_outputs_for_surface_interface_element(simulation, 4, 4, 3)
    
    def test_dirichlet_single_stage_6_plus_6(self):
        simulation = self.run_simulation('dirichlet_single_stage_6_plus_6')
        self.assert_outputs_for_surface_interface_element(simulation, 9, 3, 6)
    
    def test_dirichlet_single_stage_8_plus_8(self):
        simulation = self.run_simulation('dirichlet_single_stage_8_plus_8')
        self.assert_outputs_for_surface_interface_element(simulation, 8, 4, 3)


    #def test_neumann_multi_stage_3_plus_3(self):
    #    self.check_displacement_multi_stage('neumann_multi_stage_3_plus_3', 4)
    #
    #def test_neumann_multi_stage_4_plus_4(self):
    #    self.check_displacement_multi_stage('neumann_multi_stage_4_plus_4', 4)
    #
    #def test_neumann_multi_stage_6_plus_6(self):
    #    self.check_displacement_multi_stage('neumann_multi_stage_6_plus_6', 9)
    #
    #def test_neumann_multi_stage_8_plus_8(self):
    #    self.check_displacement_multi_stage('neumann_multi_stage_8_plus_8', 8)
    #
    #def test_dirichlet_multi_stage_3_plus_3(self):
    #    self.check_stresses_multi_stage('dirichlet_multi_stage_3_plus_3', 3)
    #
    #def test_dirichlet_multi_stage_4_plus_4(self):
    #    self.check_stresses_multi_stage('dirichlet_multi_stage_4_plus_4', 4)
    #
    #def test_dirichlet_multi_stage_6_plus_6(self):
    #    self.check_stresses_multi_stage('dirichlet_multi_stage_6_plus_6', 3)
    #
    #def test_dirichlet_multi_stage_8_plus_8(self):
    #    self.check_stresses_multi_stage('dirichlet_multi_stage_8_plus_8', 4)


    @staticmethod
    def run_simulation(test_name):
        file_path = test_helper.get_file_path(os.path.join('surface_interface_elements', test_name))
        return test_helper.run_kratos(file_path)
        
    def assert_outputs_for_surface_interface_element(self, simulation, number_of_nodes, number_of_integration_points, number_of_elements):
        displacements = test_helper.get_displacement(simulation)
        
        # Check the element at XZ
        shear_traction = -667.0
        expected_shear_displacement = shear_traction / self.shear_stiffness
        normal_traction = -333.0
        expected_normal_displacement = normal_traction / self.normal_stiffness

        ###for index in range(number_of_nodes, 2*number_of_nodes):
        ###    self.assertAlmostEqual(displacements[index][0], expected_shear_displacement)
        ###    self.assertAlmostEqual(displacements[index][1], expected_normal_displacement)
        ###    self.assertAlmostEqual(displacements[index][2], expected_shear_displacement)
        ###
        ###tractions = test_helper.get_on_integration_points(simulation, Kratos.CAUCHY_STRESS_VECTOR)
        ###tractions_horizontal_element = tractions[1]
        ###for index in range(number_of_integration_points):
        ###    self.assertAlmostEqual(tractions_horizontal_element[index][0], shear_traction)
        ###    #self.assertAlmostEqual(tractions_horizontal_element[index][1], -normal_traction)
        ###    self.assertAlmostEqual(tractions_horizontal_element[index][2], -shear_traction)
        ###    
        ###relative_displacements = test_helper.get_on_integration_points(simulation, Kratos.STRAIN)
        ###relative_displacements_horizontal_element = relative_displacements[0]
        ###for index in range(number_of_integration_points):
        ###    self.assertAlmostEqual(relative_displacements_horizontal_element[index][0], expected_shear_displacement)
        ###    self.assertAlmostEqual(relative_displacements_horizontal_element[index][1], expected_normal_displacement)
        ###    self.assertAlmostEqual(relative_displacements_horizontal_element[index][2], expected_shear_displacement)

        # Check the element at XY
        shear_traction = -667.0
        expected_shear_displacement = shear_traction / self.shear_stiffness
        normal_traction = -333.0
        expected_normal_displacement = normal_traction / self.normal_stiffness
        
        for index in range(5*number_of_nodes, 6*number_of_nodes):
            self.assertAlmostEqual(displacements[index][0], expected_shear_displacement)
            self.assertAlmostEqual(displacements[index][1], expected_shear_displacement)
            self.assertAlmostEqual(displacements[index][2], expected_normal_displacement)
        
        tractions = test_helper.get_on_integration_points(simulation, Kratos.CAUCHY_STRESS_VECTOR)
        node_index = number_of_elements * 2 // 3
        tractions_horizontal_element = tractions[node_index]
        for index in range(number_of_integration_points):
            self.assertAlmostEqual(tractions_horizontal_element[index][0], normal_traction)
            self.assertAlmostEqual(tractions_horizontal_element[index][1], shear_traction)
            self.assertAlmostEqual(tractions_horizontal_element[index][2], shear_traction)
            
        relative_displacements = test_helper.get_on_integration_points(simulation, Kratos.STRAIN)
        relative_displacements_horizontal_element = relative_displacements[node_index]
        for index in range(number_of_integration_points):
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][0], expected_normal_displacement)
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][1], expected_shear_displacement)
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][2], expected_shear_displacement)
  
        # Check the element at YZ
        shear_traction = -667.0
        expected_shear_displacement = shear_traction / self.shear_stiffness
        normal_traction = -333.0
        expected_normal_displacement = normal_traction / self.normal_stiffness
        
        for index in range(3*number_of_nodes, 4*number_of_nodes):
            print(index)
            self.assertAlmostEqual(displacements[index][0], expected_normal_displacement)
            self.assertAlmostEqual(displacements[index][1], expected_shear_displacement)
            self.assertAlmostEqual(displacements[index][2], expected_shear_displacement)
        
        #tractions = test_helper.get_on_integration_points(simulation, Kratos.CAUCHY_STRESS_VECTOR)
        #node_index = number_of_elements // 3
        #tractions_horizontal_element = tractions[node_index]
        #for index in range(number_of_integration_points):
        #    self.assertAlmostEqual(tractions_horizontal_element[index][0], normal_traction)
        #    self.assertAlmostEqual(tractions_horizontal_element[index][1], shear_traction)
        #    self.assertAlmostEqual(tractions_horizontal_element[index][2], shear_traction)
        #    
        #relative_displacements = test_helper.get_on_integration_points(simulation, Kratos.STRAIN)
        #relative_displacements_horizontal_element = relative_displacements[node_index]
        #for index in range(number_of_integration_points):
        #    self.assertAlmostEqual(relative_displacements_horizontal_element[index][0], expected_normal_displacement)
        #    self.assertAlmostEqual(relative_displacements_horizontal_element[index][1], expected_shear_displacement)
        #    self.assertAlmostEqual(relative_displacements_horizontal_element[index][2], expected_shear_displacement)
  
        # Check the element at ZX
        shear_traction = -667.0
        expected_shear_displacement = shear_traction / self.shear_stiffness
        normal_traction = -333.0
        expected_normal_displacement = normal_traction / self.normal_stiffness
        
        for index in range(number_of_nodes, 2*number_of_nodes):
            self.assertAlmostEqual(displacements[index][0], expected_shear_displacement)
            self.assertAlmostEqual(displacements[index][1], expected_normal_displacement)
            self.assertAlmostEqual(displacements[index][2], expected_shear_displacement)
        
        #tractions = test_helper.get_on_integration_points(simulation, Kratos.CAUCHY_STRESS_VECTOR)
        #node_index = 0
        #tractions_horizontal_element = tractions[node_index]
        #for index in range(number_of_integration_points):
        #    self.assertAlmostEqual(tractions_horizontal_element[index][0], normal_traction)
        #    self.assertAlmostEqual(tractions_horizontal_element[index][1], shear_traction)
        #    self.assertAlmostEqual(tractions_horizontal_element[index][2], shear_traction)
            
        #relative_displacements = test_helper.get_on_integration_points(simulation, Kratos.STRAIN)
        #relative_displacements_horizontal_element = relative_displacements[node_index]
        #for index in range(number_of_integration_points):
        #    self.assertAlmostEqual(relative_displacements_horizontal_element[index][0], expected_normal_displacement)
        #    self.assertAlmostEqual(relative_displacements_horizontal_element[index][1], expected_shear_displacement)
        #    self.assertAlmostEqual(relative_displacements_horizontal_element[index][2], expected_shear_displacement)
  
  
  
  
  
  
  
  
if __name__ == '__main__':
    KratosUnittest.main()
