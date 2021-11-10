#from KratosMultiphysics import * as Kratos
import sys
import os

sys.path.append(os.path.join('..', '..', '..'))
sys.path.append(os.path.join('..','python_scripts'))

import KratosMultiphysics as Kratos
from KratosMultiphysics.GeoMechanicsApplication import *

import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_helper

class KratosGeoMechanicsExcavationTests(KratosUnittest.TestCase):

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_base(self):

        test_name = 'base_test'
        file_path = test_helper.get_file_path(os.path.join('.','excavation_tests',test_name +'.gid'))

        simulation = test_helper.run_kratos(file_path)

        displacements = test_helper.get_displacement(simulation)

        bottom_nodes = [5, 7, 9, 12]

        # check if bottom nodes are not moving
        for bottom_node in bottom_nodes:
            self.assertEqual(0, displacements[bottom_node-1][0])
            self.assertEqual(0,displacements[bottom_node - 1][1])
            self.assertEqual(0,displacements[bottom_node - 1][2])

        top_nodes = [1, 2, 6, 10, 13]
        top_y_displacement = [3.31258e-07, 3.10315e-07, 3.01343e-07, 3.10315e-07, 3.31258e-07]
        for i in range(len(top_nodes)):
            self.assertAlmostEqual(top_y_displacement[i], displacements[top_nodes[i]-1][1])

    def test_excavation(self):
        """
        Test of a two surface geometry where 1 surface is deactivated. Values are not checked on analytical value,
        but on the fact that equal displacements should occur at both sides of the geometry.
        :return:
        """
        test_name = 'excavation_test'
        file_path = test_helper.get_file_path(os.path.join('.','excavation_tests',test_name +'.gid'))

        simulation = test_helper.run_kratos(file_path)

        displacements = test_helper.get_displacement(simulation)

        bottom_nodes = [5, 7, 9, 12]
        # check if bottom nodes are not moving
        for bottom_node in bottom_nodes:
            self.assertEqual(0, displacements[bottom_node - 1][0])
            self.assertEqual(0,displacements[bottom_node - 1][1])
            self.assertEqual(0,displacements[bottom_node - 1][2])

        top_nodes = [1,2,6,10,13]
        top_y_displacement = [0, 0, 3.16312e-07, 3.14736e-07, 3.16312e-07]
        for i in range(len(top_nodes)):
            self.assertAlmostEqual(top_y_displacement[i], displacements[top_nodes[i]-1][1])

    def test_excavation2(self):
        """
        Test of a two surface geometry where 1 surface is deactivated, on the boundary between the 2 surfaces,
        a no deformation condition is applied. Values are not checked on analytical value, but on the fact that at the
        boundary, no deformations should occur, opposite of the boundary deformations should occur.
        :return:
        """
        test_name = 'excavation_test2'
        file_path = test_helper.get_file_path(os.path.join('.','excavation_tests',test_name +'.gid'))

        simulation = test_helper.run_kratos(file_path)

        displacements = test_helper.get_displacement(simulation)
        coordinates = test_helper.get_nodal_coordinates(simulation)

        bottom_nodes = [5, 7, 9, 12]
        # check if bottom nodes are not moving
        for bottom_node in bottom_nodes:
            self.assertEqual(0, displacements[bottom_node-1][0])
            self.assertEqual(0,displacements[bottom_node - 1][1])
            self.assertEqual(0,displacements[bottom_node - 1][2])

        top_nodes = [1,2,6,10,13]
        top_y_displacement = [0, 0, 0, 2.8486711671144066e-07, 3.3932412582166974e-07]
        for i in range(len(top_nodes)):
            self.assertAlmostEqual(top_y_displacement[i], displacements[top_nodes[i]-1][1])

    def test_excavation3(self):
        """
        Test of a two surface geometry where an excavation is applied to 1 surface, but the excavation is deactivated.
        Values are not checked on analytical value, but on the fact that all nodes should have a value equal to the base test.
        :return:
        """
        test_name = 'excavation_test3'
        file_path = test_helper.get_file_path(os.path.join('.','excavation_tests',test_name +'.gid'))
        simulation = test_helper.run_kratos(file_path)

        displacements = test_helper.get_displacement(simulation)

        bottom_nodes = [5, 7, 9, 12]
        # check if bottom nodes are not moving
        for bottom_node in bottom_nodes:
            self.assertEqual(0, displacements[bottom_node-1][0])
            self.assertEqual(0,displacements[bottom_node - 1][1])
            self.assertEqual(0,displacements[bottom_node - 1][2])

        top_nodes = [1, 2, 6, 10, 13]
        top_y_displacement = [3.31258e-07, 3.10315e-07, 3.01343e-07, 3.10315e-07, 3.31258e-07]
        for i in range(len(top_nodes)):
            self.assertAlmostEqual(top_y_displacement[i], displacements[top_nodes[i]-1][1])


    @KratosUnittest.skip("unit test skipped as it is not ready")
    def test_phases(self):
        """
        Test of a two surface geometry where 1 surface is deactivated. Values are not checked on analytical value,
        but on the fact that equal displacements should occur at both sides of the geometry.
        :return:
        """
        print("todo")
        # test_name_1 = 'excavation_test'
        # test_name_2 = 'excavation_test4'
        # file_path_1 = GetFilePath(os.path.join('test_data',test_name_1))
        # file_path_2 = GetFilePath(os.path.join('test_data', test_name_2))
        #
        # parameter_file_name_1 = os.path.join(file_path_1, 'ProjectParameters.json')
        # parameter_file_name_2 = os.path.join(file_path_2, 'ProjectParameters.json')
        #
        # os.chdir(file_path_1)
        # with open(parameter_file_name_1, 'r') as parameter_file:
        #     parameters_1 = Kratos.Parameters(parameter_file.read())
        #
        # os.chdir(file_path_2)
        # with open(parameter_file_name_2, 'r') as parameter_file:
        #     parameters_2 = Kratos.Parameters(parameter_file.read())
        #
        # model = Kratos.Model()
        #
        # stage_1 = analysis.GeoMechanicsAnalysis(model, parameters_1)
        # stage_2 = analysis.GeoMechanicsAnalysis(model, parameters_2)
        #
        # stage_1.Run()
        # stage_2.Run()
        #
        # displacements = get_displacement(stage_2)
        #
        # bottom_nodes = [5, 7, 9, 12]
        # # check if bottom nodes are not moving
        # for bottom_node in bottom_nodes:
        #     self.assertEqual(0, displacements[bottom_node - 1][0])
        #     self.assertEqual(0,displacements[bottom_node - 1][1])
        #     self.assertEqual(0,displacements[bottom_node - 1][2])
        #
        # top_nodes = [1, 2, 6, 10, 13]
        # top_y_displacement = [0, 0, 3.16312e-07, 3.14736e-07, 3.16312e-07]
        # for i in range(len(top_nodes)):
        #     self.assertAlmostEqual(top_y_displacement[i], displacements[top_nodes[i]-1][1])

    def testNightlyFirstExample(self):
        self.assertEqual(True, True)

    def testNightlySecondExample(self):
        self.assertEqual(True, True)

if __name__ == '__main__':
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([KratosGeoMechanicsExcavationTests]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)
