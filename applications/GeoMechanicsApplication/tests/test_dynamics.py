import sys
import os

sys.path.append(os.path.join('..', '..', '..'))
sys.path.append(os.path.join('..', 'python_scripts'))

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsDynamicsTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with other FE software
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    @KratosUnittest.skip("This test is very long and should be shortend")
    def test_wave_through_drained_linear_elastic_soil(self):
        """
        Test dynamic calculation on a drained linear elastic soil column. a line load of -1kN is instantly placed
        on the soil column. The soil parameters are chosen such that after 0.002 seconds, the wave is reflected at the
        bottom of the geometry such that half the stress in the soil column is cancelled out.
        :return:
        """
        test_name = 'test_1d_wave_prop_drained_soil.gid'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))

        simulation = test_helper.run_kratos(file_path)

        # get effective stress
        efective_stresses = test_helper.get_cauchy_stress_tensor(simulation)
        efective_stresses_yy = [integration_point[1,1] for element in efective_stresses for integration_point in element]

        # get coordinates of the gauss points
        gauss_coordinates = test_helper.get_gauss_coordinates(simulation)
        gauss_coordinates_y = [integration_point[1] for element in gauss_coordinates for integration_point in element]

        # calculate the expected effective stress in the soil column
        expected_effective_stresses_yy = [-1000 if gauss_coordinate_y >= -0.49 else 0 if gauss_coordinate_y <= -0.51
                                          else -500 for gauss_coordinate_y in gauss_coordinates_y]

        # calculate root mean square error
        square_errors = [(efective_stress_yy - expected_effective_stresses_yy[idx])**2
                         for idx, efective_stress_yy in enumerate(efective_stresses_yy)]
        rmse = (sum(square_errors)/len(square_errors))**0.5

        # assert root mean square error, the allowable error is 10% of the applied load
        self.assertLess(rmse, 100)


    @KratosUnittest.skip("unit test skipped as it is not ready")
    def test_wave_through_undrained_linear_elastic_soil(self):
        test_name = 'test_1d_confined_undrained_wave.gid'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        simulation = test_helper.run_kratos(file_path)
        pass

if __name__ == '__main__':
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([KratosGeoMechanicsDynamicsTests]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)
