import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis


class KratosGeoMechanicsMohrCoulombElastoPlasticTangentMatrixTests(
        KratosUnittest.TestCase):

    def setUp(self):
        self.work_folder = os.path.join(
            "test_mohr_coulomb_elastoplastic_tangent_matrix", "common")
        self.print_reference_values = True

    def test_triangle_under_gravity(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            with open("ProjectParameters.json", "r") as parameter_file:
                parameters = Kratos.Parameters(parameter_file.read())

            if self.print_reference_values:
                self._add_reference_values_output(parameters)
            else:
                self._add_reference_values_check(parameters)

            model = Kratos.Model()
            simulation = analysis.GeoMechanicsAnalysis(model, parameters)
            simulation.Run()

    @staticmethod
    def _add_reference_values_output(parameters):
        json_output_settings = Kratos.Parameters("""{
            "python_module": "json_output_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "JsonOutputProcess",
            "Parameters": {
                "model_part_name": "PorousDomain",
                "output_file_name": "triangle_under_gravity_results.json",
                "output_variables": ["DISPLACEMENT", "WATER_PRESSURE"],
                "time_frequency": 1.0
            }
        }""")
        parameters["processes"].AddEmptyArray("json_check_process_list")
        parameters["processes"]["json_check_process_list"].Append(
            json_output_settings)

    @staticmethod
    def _add_reference_values_check(parameters):
        json_check_settings = Kratos.Parameters("""{
            "python_module": "from_json_check_result_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "FromJsonCheckResultProcess",
            "Parameters": {
                "model_part_name": "PorousDomain",
                "input_file_name": "triangle_under_gravity_results.json",
                "check_variables": ["DISPLACEMENT", "WATER_PRESSURE"],
                "tolerance": 1.0e-7,
                "relative_tolerance": 1.0e-5,
                "time_frequency": 1.0
            }
        }""")
        parameters["processes"].AddEmptyArray("json_check_process_list")
        parameters["processes"]["json_check_process_list"].Append(
            json_check_settings)


if __name__ == "__main__":
    KratosUnittest.main()
