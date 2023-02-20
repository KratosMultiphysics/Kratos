import os
import sys
import shutil

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsParameterFieldTests(KratosUnittest.TestCase):
    """
    This class contains tests which check if the lysmer absorbing boundary works in a 1d column made out of different
    element types
    """

    # def setUp(self):
    #     self.E = 10000
    #     self.nu = 0.2
    #     self.rho = 2.65 * 0.7
    #     self.load = -10
    #     self.height_column = 10
    #     self.vp = None
    #     self.expected_velocity = None
    #
    #     self.calculate_expected_p_wave_velocity()
    #     self.calculate_expected_velocity_column()
    #
    # def tearDown(self):
    #     # Code here will be placed AFTER every test in this TestCase.
    #     pass

    def test_parameter_field_with_function(self):
        """
        Tests the lysmer absorbing boundary condition on a column made of rectangulars. The boundary is a 2d2n line.

        :return:
        """
        test_name = os.path.join("test_parameter_field", "parameter_field_input")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # get element centers
        elements = simulation._list_of_output_processes[0].model_part.Elements
        center_coords = [element.GetGeometry().Center() for element in elements]

        results = test_helper.get_on_integration_points(simulation, Kratos.YOUNG_MODULUS)

        for center_coord, res in zip(center_coords, results):
            expected_res = 20000 * center_coord[0] + 30000 * center_coord[1]
            self.assertAlmostEqual(expected_res, res[0])


    def test_parameter_field_with_python(self):
        """
        Tests the lysmer absorbing boundary condition on a column made of rectangulars. The boundary is a 2d2n line.

        :return:
        """
        test_name = os.path.join("test_parameter_field", "parameter_field_python")
        file_path = test_helper.get_file_path(test_name)

        custom_python_file = os.path.join(file_path, "custom_field.py")

        new_custom_script_path = os.path.join(os.path.dirname(KratosGeo.__file__), "user_defined_scripts")
        shutil.copy(custom_python_file,new_custom_script_path)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # get element centers
        elements = simulation._list_of_output_processes[0].model_part.Elements
        center_coords = [element.GetGeometry().Center() for element in elements]

        results = test_helper.get_on_integration_points(simulation, Kratos.YOUNG_MODULUS)

        for center_coord, res in zip(center_coords, results):
            expected_res = 20000 * center_coord[0] + 30000 * center_coord[1]
            self.assertAlmostEqual(expected_res, res[0])

    def test_parameter_field_with_json(self):
        """
        Tests the lysmer absorbing boundary condition on a column made of rectangulars. The boundary is a 2d2n line.

        :return:
        """
        test_name = os.path.join("test_parameter_field", "parameter_field_json")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # get element centers
        elements = simulation._list_of_output_processes[0].model_part.Elements
        center_coords = [element.GetGeometry().Center() for element in elements]

        results = test_helper.get_on_integration_points(simulation, Kratos.YOUNG_MODULUS)

        for center_coord, res in zip(center_coords, results):
            expected_res = 20000 * center_coord[0] + 30000 * center_coord[1]
            self.assertAlmostEqual(expected_res, res[0])



if __name__ == '__main__':
    KratosUnittest.main()
