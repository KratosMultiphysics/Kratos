import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsPipingMethodTests(KratosUnittest.TestCase):

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_self_recursive_solution_step(self):
        test_name = 'test_piping_method'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        simulation = test_helper.run_kratos(file_path)
        model_part = simulation._list_of_output_processes[0].model_part
        self.assertTrue(model_part is not None)


if __name__ == '__main__':
    suites = KratosUnittest.KratosSuites
    small_suite = suites['small'] # These tests are executed by the continuous integration tool
    small_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([KratosGeoMechanicsPipingMethodTests]))
    all_suite = suites['all']
    all_suite.addTests(small_suite)
    KratosUnittest.runTests(suites)
