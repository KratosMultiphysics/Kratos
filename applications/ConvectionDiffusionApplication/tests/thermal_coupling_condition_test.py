import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtilities
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis

class ConvectionDiffusionAnalysisThermalCouplingCondition(ConvectionDiffusionAnalysis):

    def __init__(self, model, project_parameters, transfer_coefficient):
        super().__init__(model,project_parameters)
        self._transfer_coefficient = transfer_coefficient

    def ModifyAfterSolverInitialize(self):
        main_model_part = self.model.GetModelPart("ThermalModelPart")
        interface_prop = main_model_part.GetProperties(2,0)
        interface_prop.SetValue(KratosMultiphysics.ConvectionDiffusionApplication.TRANSFER_COEFFICIENT, self._transfer_coefficient)

class ThermalCouplingConditionTest(UnitTest.TestCase):

    def setUp(self):
        self.work_folder = "thermal_coupling_condition_test"
        self.check_absolute_tolerance = 1.0e-7
        self.check_relative_tolerance = 1.0e-5
        self.print_output = False
        self.print_reference_values = False

    def runTest(self):
        with UnitTest.WorkFolderScope(self.work_folder,__file__):
            with open("ProjectParameters.json",'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())
                model = KratosMultiphysics.Model()

                if self.print_output:
                    self._AddOutput()

                if self.print_reference_values:
                    self._AddReferenceValuesOutput()
                else:
                    self._AddReferenceValuesCheck()

                # running
                self.simulation = ConvectionDiffusionAnalysisThermalCouplingCondition(model, self.parameters, self._transfer_coefficient)
                self.simulation.Run()

    def tearDown(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            KratosUtilities.DeleteFileIfExisting('.time')

    def testThermalCouplingCondition(self):
        self._transfer_coefficient = 1.0
        self.setUp()
        self.runTest()
        self.tearDown()

    def _AddOutput(self):
        gid_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "help"          : "This process writes postprocessing files for GiD",
            "Parameters"    : {
                "model_part_name"        : "ThermalModelPart",
                "output_name"            : "TO_BE_DEFINED",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"               : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "SingleFile"
                        },
                        "file_label"                  : "step",
                        "output_control_type"         : "step",
                        "output_interval"             : 1.0,
                        "body_output"                 : true,
                        "node_output"                 : false,
                        "skin_output"                 : false,
                        "plane_output"                : [],
                        "nodal_results"               : ["TEMPERATURE"]
                    },
                    "point_data_configuration"  : []
                }
            }
        }""")
        output_name = "thermal_coupling_condition_test"
        gid_output_settings["Parameters"]["output_name"].SetString(output_name)
        self.parameters["output_processes"]["gid_output"].Append(gid_output_settings)

    def _AddReferenceValuesOutput(self):
        json_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "json_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "JsonOutputProcess",
            "Parameters"    : {
                "output_variables" : ["TEMPERATURE"],
                "output_file_name" : "TO_BE_DEFINED",
                "model_part_name"  : "ThermalModelPart",
                "time_frequency"   : 1.0
            }
        }""")
        output_file_name = "thermal_coupling_condition_test_results.json"
        json_output_settings["Parameters"]["output_file_name"].SetString(output_file_name)
        self.parameters["processes"]["json_check_process_list"].Append(json_output_settings)

    def _AddReferenceValuesCheck(self):
        json_check_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "from_json_check_result_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "FromJsonCheckResultProcess",
            "Parameters"    : {
                "check_variables"      : ["TEMPERATURE"],
                "input_file_name"      : "TO_BE_DEFINED",
                "model_part_name"      : "ThermalModelPart",
                "tolerance"            : 0.0,
                "relative_tolerance"   : 0.0,
                "time_frequency"       : 1.0
            }
        }""")
        input_file_name = "thermal_coupling_condition_test_results.json"
        json_check_settings["Parameters"]["input_file_name"].SetString(input_file_name)
        json_check_settings["Parameters"]["tolerance"].SetDouble(self.check_absolute_tolerance)
        json_check_settings["Parameters"]["relative_tolerance"].SetDouble(self.check_relative_tolerance)
        self.parameters["processes"]["json_check_process_list"].Append(json_check_settings)


if __name__ == "__main__":
    UnitTest.main()
