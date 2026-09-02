import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis_rve import FluidDynamicsAnalysisRve

class TestFluidRVETest(KratosUnittest.TestCase):

    def setUp(self):
        # Change to True to generate GiD Output
        self.print_output = False

        # Change to True if no reference value file exists
        self.check_tolerance = 1e-6
        self.check_relative_tolerance = 1e-8
        self.print_reference_values = False

    def test_fluid_rve_computation_2d(self):
        #Within location context:
        with KratosUnittest.WorkFolderScope(".",__file__):
            with open("FluidRVETest/fluid_rve_test_parameters.json", 'r') as parameter_file:
                self.parameters =  KratosMultiphysics.Parameters(parameter_file.read())

            self.parameters["solver_settings"]["domain_size"].SetInt(2)
            self.parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString("FluidRVETest/fluid_rve_test_2D")
            self.parameters["solver_settings"]["material_import_settings"]["materials_filename"].SetString("FluidRVETest/fluid_rve_test_materials_2D.json")
            self.parameters["rve_settings"]["boundary_mp_name"].SetString("FluidModelPart.Slip2D.Boundaries")
            self.parameters["solver_settings"]["skin_parts"][0].SetString("Slip2D")
            self.parameters["processes"]["initial_conditions_process_list"][0]["Parameters"]["value"][0].SetInt(0)
            self.parameters["processes"]["initial_conditions_process_list"][0]["Parameters"]["value"][1].SetInt(0)

            self._aux_fluid_rve_computation()

    def test_fluid_rve_computation_3d(self):
        #Within location context:
        with KratosUnittest.WorkFolderScope(".",__file__):
            with open("FluidRVETest/fluid_rve_test_parameters.json", 'r') as parameter_file:
                self.parameters =  KratosMultiphysics.Parameters(parameter_file.read())

            self.parameters["solver_settings"]["domain_size"].SetInt(3)
            self.parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString("FluidRVETest/fluid_rve_test_3D")
            self.parameters["solver_settings"]["material_import_settings"]["materials_filename"].SetString("FluidRVETest/fluid_rve_test_materials_3D.json")
            self.parameters["rve_settings"]["boundary_mp_name"].SetString("FluidModelPart.Slip3D.Boundaries")
            self.parameters["solver_settings"]["skin_parts"][0].SetString("Slip3D")
            self.parameters["processes"]["initial_conditions_process_list"][0]["Parameters"]["value"][0].SetInt(0)
            self.parameters["processes"]["initial_conditions_process_list"][0]["Parameters"]["value"][1].SetInt(0)
            self.parameters["processes"]["initial_conditions_process_list"][0]["Parameters"]["value"][2].SetInt(0)

            self._aux_fluid_rve_computation()

    def _aux_fluid_rve_computation(self):
        # If required, add the output process to the test settings
        if self.print_output:
            self._AddOutput()

        # If required, add the reference values output process to the test settings
        if self.print_reference_values:
            self._AddReferenceValuesOutput()
        else:
            self._AddReferenceValuesCheck()

        # Run RVE simulation
        model = KratosMultiphysics.Model()
        self.simulation = FluidDynamicsAnalysisRve(model, self.parameters)
        self.simulation.Run()

    def _AddOutput(self):
        domain_size = self.parameters["solver_settings"]["domain_size"].GetInt()
        output_settings = KratosMultiphysics.Parameters(R'''[{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "help"          : "This process writes postprocessing files for GiD",
            "Parameters"    : {
                "model_part_name"        : "FluidModelPart",
                "output_name"            : "TO_BE_SET",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"         : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "SingleFile"
                        },
                        "file_label"          : "time",
                        "output_control_type" : "step",
                        "output_interval"     : 1,
                        "body_output"         : true,
                        "node_output"         : false,
                        "skin_output"         : false,
                        "plane_output"        : [],
                        "nodal_results"       : ["VELOCITY","PRESSURE"],
                        "gauss_point_results" : [],
                        "nodal_flags_results": ["MASTER","SLAVE"]
                    },
                    "point_data_configuration"  : []
                }
            }
        }]''')
        output_settings[0]["Parameters"]["output_name"].SetString("FluidRVETest/fluid_rve_test_{}D".format(domain_size))
        self.parameters["output_processes"].AddValue("gid_output", output_settings)

    def _AddReferenceValuesOutput(self):
        json_output_settings = KratosMultiphysics.Parameters("""{
            "name" : "Processes.KratosMultiphysics.JsonOutputProcess",
            "Parameters"    : {
                "output_variables" : ["VELOCITY","PRESSURE"],
                "output_file_name" : "FluidRVETest/fluid_rve_test_results",
                "model_part_name"  : "FluidModelPart",
                "time_frequency"   : 0.1
            }
        }""")
        domain_size = self.parameters["solver_settings"]["domain_size"].GetInt()
        output_file_name = json_output_settings["Parameters"]["output_file_name"].GetString()
        output_file_name += "_{0}D.json".format(domain_size)
        json_output_settings["Parameters"]["output_file_name"].SetString(output_file_name)
        self.parameters["processes"]["json_check_process_list"].Append(json_output_settings)

    def _AddReferenceValuesCheck(self):
        json_check_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "from_json_check_result_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "FromJsonCheckResultProcess",
            "Parameters"    : {
                "check_variables"      : ["VELOCITY_X","VELOCITY_Y","PRESSURE"],
                "input_file_name"      : "FluidRVETest/fluid_rve_test_results",
                "model_part_name"      : "FluidModelPart",
                "tolerance"            : 0.0,
                "relative_tolerance"   : 0.0,
                "time_frequency"       : 0.1
            }
        }""")
        domain_size = self.parameters["solver_settings"]["domain_size"].GetInt()
        input_file_name = json_check_settings["Parameters"]["input_file_name"].GetString()
        input_file_name += "_{0}D.json".format(domain_size)
        json_check_settings["Parameters"]["input_file_name"].SetString(input_file_name)
        json_check_settings["Parameters"]["tolerance"].SetDouble(self.check_tolerance)
        json_check_settings["Parameters"]["relative_tolerance"].SetDouble(self.check_relative_tolerance)
        self.parameters["processes"]["json_check_process_list"].Append(json_check_settings)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

