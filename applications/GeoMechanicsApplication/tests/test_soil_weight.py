import sys
import os

sys.path.append(os.path.join('..', '..', '..'))
sys.path.append(os.path.join('..', 'python_scripts'))

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsSoilWeightTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_saturated_below_phreatic_quad_8n(self):
        test_name = 'test_soil_weight_saturated_below_phreatic_quad_8n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        self.assert_total_stress(simulation, -32.886751345948134)

    def test_saturated_below_phreatic_quad_4n(self):
        test_name = 'test_soil_weight_saturated_below_phreatic_quad_4n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        self.assert_total_stress(simulation, -31.44337567297404)

    def test_saturated_below_and_above_phreatic_quad_8n(self):
        test_name = 'test_soil_weight_saturated_below_and_above_phreatic_quad_8n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        self.assert_total_stress(simulation, -37.886751345948134)

    def test_saturated_below_and_above_phreatic_quad_4n(self):
        test_name = 'test_soil_weight_saturated_below_and_above_phreatic_quad_4n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        self.assert_total_stress(simulation, -36.44337567297406)

    def test_van_genuchten_above_phreatic_quad_8n(self):
        test_name = 'test_soil_weight_van_genuchten_above_phreatic_quad_8n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        self.assert_total_stress(simulation, -36.65995752451169)

    def test_van_genuchten_above_phreatic_quad_4n(self):
        test_name = 'test_soil_weight_van_genuchten_above_phreatic_quad_4n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        self.assert_total_stress(simulation, -35.21658185153758)

    def assert_total_stress(self, simulation, expected_value):
        """
        Assert results of a linear elastic block. The sides of the block can move freely in vertical direction and are
        fixed in horizontal direction. The bottom of the block is fixed. A gravitational acceleration of 10 m/s2 is applied to the 
        whole body. Results are: total stresses.

        :param simulation: Kratos simulation
        :param expected_value: expected total stress at the bottom of the model.
        :return:
        """
        total_stresses = test_helper.get_total_stress_tensor(simulation)
        total_stresses_yy = [integration_point[1,1] for element in total_stresses for integration_point in element]

        # Assert integration point information
        min_total_stresses_yy = 1e10
        for idx, total_stress_xx in enumerate(total_stresses_yy):
            min_total_stresses_yy = min(min_total_stresses_yy, total_stresses_yy[idx])
            
        self.assertAlmostEqual(expected_value, min_total_stresses_yy)


if __name__ == '__main__':
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([KratosGeoMechanicsSoilWeightTests]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)
