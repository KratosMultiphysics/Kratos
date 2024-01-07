from math import sqrt

import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS

import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.RANSApplication.response_functions.utilities import GetResponseFunctionOutputProcess



class ResponseFunctionUtilitiesTest(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("FluidModelPart")
        cls.model_part.SetBufferSize(2)

        cls.model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 2)
        cls.model_part.ProcessInfo.SetValue(Kratos.STEP, 1)
        cls.model_part.ProcessInfo.SetValue(Kratos.DELTA_TIME, 0.1)

    def setUp(self):
        pass

    def testGetResponseFunctionOutputProcess(self):
        optimization_response_parameters = Kratos.Parameters("""{
            "custom_settings": {
                "magnitude_square_to_power": 1,
                "model_part_name": "GENERIC_Investigation_WindComfort_AreaFront",
                "variable_name": "VELOCITY"
            },
            "response_type": "domain_integrated_3d_vector_magnitude_square_power_mean"
        }""")

        primal_parameters = Kratos.Parameters("""{
            "output_processes": {
                "hdf5_output": [
                    {
                        "Parameters": {
                            "file_settings": {
                                "echo_level": 0,
                                "file_access_mode": "truncate",
                                "file_name": "results/primal/<model_part_name>-final.h5",
                                "max_files_to_keep": 10,
                                "time_format": "0.6f"
                            },
                            "model_part_name": "FluidModelPart",
                            "model_part_output_settings": {}
                        },
                        "kratos_module": "KratosMultiphysics.RANSApplication",
                        "python_module": "primal_hdf5_output_process"
                    }
                ],
                "response_function": [
                    {
                        "Parameters": {
                            "model_part_name": "FluidModelPart",
                            "output_file_settings": {
                                "file_name": "domain_integrated_3d_vector_magnitude_square_1_mean",
                                "output_path": "results/ascii",
                                "write_buffer_size": -1
                            },
                            "response_settings": {
                                "dof_position": 0,
                                "entities_to_consider": [
                                    "conditions"
                                ],
                                "flag_to_be_used": "STRUCTURE",
                                "magnitude_square_to_power": 3,
                                "model_part_name": "GENERIC_Investigation_WindComfort_AreaFront",
                                "start_time": 0.0,
                                "variable_name": "VELOCITY"
                            },
                            "response_type": "domain_integrated_3d_vector_magnitude_square_power_mean"
                        },
                        "kratos_module": "KratosMultiphysics.FluidDynamicsApplication",
                        "python_module": "response_function_output_process"
                    },
                    {
                        "Parameters": {
                            "model_part_name": "FluidModelPart",
                            "output_file_settings": {
                                "file_name": "domain_integrated_3d_vector_magnitude_square_1_mean",
                                "output_path": "results/ascii",
                                "write_buffer_size": -1
                            },
                            "response_settings": {
                                "dof_position": 0,
                                "entities_to_consider": [
                                    "conditions"
                                ],
                                "flag_to_be_used": "STRUCTURE",
                                "magnitude_square_to_power": 1,
                                "model_part_name": "GENERIC_Investigation_WindComfort_AreaFront",
                                "start_time": 0.0,
                                "variable_name": "VELOCITY"
                            },
                            "response_type": "domain_integrated_3d_vector_magnitude_square_power_mean"
                        },
                        "kratos_module": "KratosMultiphysics.FluidDynamicsApplication",
                        "python_module": "response_function_output_process"
                    }
                ]
            }
        }""")

        result = GetResponseFunctionOutputProcess(primal_parameters, "FluidModelPart.GENERIC_Investigation_WindComfort_AreaFront", optimization_response_parameters)
        self.assertTrue(result.IsEquivalentTo(primal_parameters["output_processes"]["response_function"][1]))

if __name__ == '__main__':
    UnitTest.main()
