import sys
import os

sys.path.append(os.path.join('..', '..', '..'))
sys.path.append(os.path.join('..', 'python_scripts'))

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsSteadyStateGroundWaterFlowTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_SteadyState_OneDFlow_2D3N(self):
        test_name = 'test_SteadyState_OneDFlow_2D3N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        self.assert_outflow_discharge(simulation, 1.0)

    def test_SteadyState_DamConfinedFlow_2D3N(self):
        test_name = 'test_SteadyState_DamConfinedFlow_2D3N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        self.assert_outflow_discharge(simulation, 1.11309991)

    def test_SteadyState_DamConfinedFlow_2D6N(self):
        test_name = 'test_SteadyState_DamConfinedFlow_2D6N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        self.assert_outflow_discharge(simulation, 1.08514372)

    def test_SteadyState_DamConfinedFlowWithSheetPile_2D6N(self):
        test_name = 'test_SteadyState_DamConfinedFlowWithSheetPile_2D6N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        self.assert_outflow_discharge(simulation, 0.827800245)


    def calculate_outflow_discharge(self, simulation):
        """
        calculate total outflow discharge.

        :param simulation: Kratos simulation
        :return: total outflow discharge
        """

        hydraulic_discharge = test_helper.get_hydraulic_discharge(simulation)

        from functools import reduce
        return reduce(lambda a,b: a+b if b>0 else a, hydraulic_discharge, 0)

    def assert_outflow_discharge(self, simulation, expected_value):
        """
        assert total outflow discharge.

        :param simulation: Kratos simulation
        :param expected_value: expected value
        :return:
        """

        outflow_discharge = self.calculate_outflow_discharge(simulation)
        error = abs(outflow_discharge - expected_value) / (expected_value + 1e-60)
        tolerated_error = 1e-6
        self.assertLess(error, tolerated_error)


if __name__ == '__main__':
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([KratosGeoMechanicsSteadyStateGroundWaterFlowTests]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)
