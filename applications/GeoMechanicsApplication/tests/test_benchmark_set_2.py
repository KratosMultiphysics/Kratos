
import sys
import os

sys.path.append(os.path.join('..', '..', '..'))
sys.path.append(os.path.join('..', 'python_scripts'))

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsBenchmarkSet2(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with other FE software
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_benchmark2_1(self):
        """
        In this benchmark a two stage model is tested, where a simple clay dike is put on top of a sand layer.
        The test is calculated using gravity loading. And UDSM-Mohr-Coulomb

        This test compares minimum, maximum stress in the last stage.

        :return:
        """
        test_name = os.path.join('Simple_Dike_Gravity_Loading','line_body_all_stage_new_units_kPa')
        project_path = test_helper.get_file_path(os.path.join('.', test_name))
        n_stages = 2
        stages = test_helper.run_stages(project_path, n_stages)

        max_x_total_stress_plaxis = 0.0
        min_x_total_stress_plaxis = -140430.9

        max_y_total_stress_plaxis = -1377.6
        min_y_total_stress_plaxis = -424253.415

        cauchy_stresses_kratos = test_helper.get_cauchy_stress_tensor(stages[-1])

        x_cauchy_stress_kratos = [[gauss_point[0,0] for gauss_point in element] for element in cauchy_stresses_kratos]
        y_cauchy_stress_kratos = [[gauss_point[1,1] for gauss_point in element] for element in cauchy_stresses_kratos]

        min_x_cauchy_stress_kratos,max_x_cauchy_stress_kratos = min(min(x_cauchy_stress_kratos)), max(max(x_cauchy_stress_kratos))
        min_y_cauchy_stress_kratos, max_y_cauchy_stress_kratos = min(min(y_cauchy_stress_kratos)), max(max(y_cauchy_stress_kratos))

        relative_precision_stress = 0.05

        absolute_precision_stress = 1e3

        self.assertAlmostEqual(max_x_total_stress_plaxis, max_x_cauchy_stress_kratos, delta=max(abs(max_x_total_stress_plaxis*relative_precision_stress),absolute_precision_stress))
        self.assertAlmostEqual(min_x_total_stress_plaxis, min_x_cauchy_stress_kratos, delta=max(abs(min_x_total_stress_plaxis*relative_precision_stress),absolute_precision_stress))

        self.assertAlmostEqual(max_y_total_stress_plaxis, max_y_cauchy_stress_kratos,  delta=max(max_y_total_stress_plaxis*relative_precision_stress,absolute_precision_stress))
        self.assertAlmostEqual(min_y_total_stress_plaxis, min_y_cauchy_stress_kratos,  delta=max(abs(min_x_total_stress_plaxis*relative_precision_stress),absolute_precision_stress))

    def test_benchmark2_2(self):
        """
        In this benchmark a two stage model is tested, where a simple clay dike is put on top of a sand layer.
        The test is calculated using gravity loading. And UMAT-Mohr-Coulomb

        This test compares minimum, maximum stress in the last stage.

        :return:
        """
        test_name = os.path.join('Simple_Dike_Gravity_Loading','simple_dike_test_with_gravity_umat.gid')
        project_path = test_helper.get_file_path(os.path.join('.', test_name))
        n_stages = 2
        stages = test_helper.run_stages(project_path, n_stages)

        max_x_total_stress_plaxis = 0.0
        min_x_total_stress_plaxis = -140430.9

        max_y_total_stress_plaxis = -1377.6
        min_y_total_stress_plaxis = -424253.415

        cauchy_stresses_kratos = test_helper.get_cauchy_stress_tensor(stages[-1])

        x_cauchy_stress_kratos = [[gauss_point[0,0] for gauss_point in element] for element in cauchy_stresses_kratos]
        y_cauchy_stress_kratos = [[gauss_point[1,1] for gauss_point in element] for element in cauchy_stresses_kratos]

        min_x_cauchy_stress_kratos,max_x_cauchy_stress_kratos = min(min(x_cauchy_stress_kratos)), max(max(x_cauchy_stress_kratos))
        min_y_cauchy_stress_kratos, max_y_cauchy_stress_kratos = min(min(y_cauchy_stress_kratos)), max(max(y_cauchy_stress_kratos))

        relative_precision_stress = 0.05

        absolute_precision_stress = 1e3

        self.assertAlmostEqual(max_x_total_stress_plaxis, max_x_cauchy_stress_kratos, delta=max(abs(max_x_total_stress_plaxis*relative_precision_stress),absolute_precision_stress))
        self.assertAlmostEqual(min_x_total_stress_plaxis, min_x_cauchy_stress_kratos, delta=max(abs(min_x_total_stress_plaxis*relative_precision_stress),absolute_precision_stress))

        self.assertAlmostEqual(max_y_total_stress_plaxis, max_y_cauchy_stress_kratos,  delta=max(max_y_total_stress_plaxis*relative_precision_stress,absolute_precision_stress))
        self.assertAlmostEqual(min_y_total_stress_plaxis, min_y_cauchy_stress_kratos,  delta=max(abs(min_x_total_stress_plaxis*relative_precision_stress),absolute_precision_stress))

    @KratosUnittest.skip("unit test skipped as it is not ready")
    def test_benchmark2_3(self):
        """
        In this benchmark a four stage model is tested.
        Stage 1: initialisation sand layer
        Stage 2: addition of sheetpile
        Stage 3: addition of dike
        Stage 4: addition of water

        The test is calculated using gravity loading.

        This test compares minimum, maximum stress and displacement in the third and fourth stage.

        :return:
        """

        test_name = r'dike_with_sheetpile_all_stage'
        project_path = test_helper.get_file_path(os.path.join('.', test_name))
        n_stages = 3
        stages = test_helper.run_stages(project_path, n_stages)

        # max_x_disp_plaxis = 0.01279
        # min_x_disp_plaxis = -0.012234
        #
        # max_y_disp_plaxis = 0.0
        # min_y_disp_plaxis = -0.21458
        #
        # max_x_total_stress_plaxis = -311
        # min_x_total_stress_plaxis = -146256
        #
        # max_y_total_stress_plaxis = -674
        # min_y_total_stress_plaxis = -432068

if __name__ == '__main__':
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([KratosGeoMechanicsBenchmarkSet2]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)
