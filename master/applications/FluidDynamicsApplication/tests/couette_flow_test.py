# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

class CouetteFlowTest(KratosUnittest.TestCase):
    def testCouetteFlow2DSymbolicStokes(self):
        self.solver_type = "MonolithicStokes"
        self.element_type = "symbolic_stokes"
        self.time_scheme = "bdf2"
        self._CustomizeTestSettings()
        self._RunCouetteFlowTest()

    def testCouetteFlow2DWeaklyCompressibleNavierStokes(self):
        self.solver_type = "Monolithic"
        self.element_type = "weakly_compressible"
        self.time_scheme = "bdf2"
        self._CustomizeTestSettings()
        self._RunCouetteFlowTest()

    def setUp(self):
        self.print_output = False
        self.check_tolerance = 1.0e-6
        self.check_relative_tolerance = 1.0e-8
        self.print_reference_values = False
        self.work_folder = "CouetteFlowTest"
        settings_filename = "ProjectParameters.json"

        # Read the simulation settings
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            with open(settings_filename,'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())

    def _RunCouetteFlowTest(self):
        # Create the test simulation
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            self.model = KratosMultiphysics.Model()
            simulation = FluidDynamicsAnalysis(self.model, self.parameters)
            simulation.Run()

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            KratosUtilities.DeleteFileIfExisting('couette_flow_test.time')

    def _CustomizeTestSettings(self):
        # Customize simulation settings
        self.parameters["solver_settings"]["solver_type"].SetString(self.solver_type)
        self.parameters["solver_settings"]["formulation"]["element_type"].SetString(self.element_type)
        self.parameters["solver_settings"]["time_scheme"].SetString(self.time_scheme)

        # If required, add the output process to the test settings
        if self.print_output:
            self._AddOutput()

        # If required, add the reference values output process to the test settings
        if self.print_reference_values:
            self._AddReferenceValuesOutput()
        else:
            self._AddReferenceValuesCheck()

    def _AddOutput(self):
        gid_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "help"          : "This process writes postprocessing files for GiD",
            "Parameters"    : {
                "model_part_name"        : "FluidModelPart.fluid_computational_model_part",
                "output_name"            : "couette_flow_test",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"               : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "SingleFile"
                        },
                        "file_label"                  : "time",
                        "output_control_type"         : "step",
                        "output_interval"             : 1,
                        "body_output"                 : true,
                        "node_output"                 : false,
                        "skin_output"                 : false,
                        "plane_output"                : [],
                        "nodal_results"               : ["VELOCITY","PRESSURE"],
                        "gauss_point_results"         : [],
                        "nodal_nonhistorical_results" : []
                    },
                    "point_data_configuration"  : []
                }
            }
        }""")
        output_name = gid_output_settings["Parameters"]["output_name"].GetString()
        output_name += "_" + self.element_type
        gid_output_settings["Parameters"]["output_name"].SetString(output_name)
        self.parameters["output_processes"]["gid_output"].Append(gid_output_settings)

    def _AddReferenceValuesOutput(self):
        json_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "json_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "JsonOutputProcess",
            "Parameters"    : {
                "output_variables" : ["VELOCITY","PRESSURE"],
                "output_file_name" : "reference_couette_flow_test",
                "model_part_name"  : "FluidModelPart.FluidParts_Fluid",
                "time_frequency"   : 99.0
            }
        }""")
        output_file_name = json_output_settings["Parameters"]["output_file_name"].GetString()
        output_file_name += "_" + self.element_type + ".json"
        json_output_settings["Parameters"]["output_file_name"].SetString(output_file_name)
        self.parameters["processes"]["json_check_process_list"].Append(json_output_settings)

    def _AddReferenceValuesCheck(self):
        json_check_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "from_json_check_result_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "FromJsonCheckResultProcess",
            "Parameters"    : {
                "check_variables"      : ["VELOCITY_X","VELOCITY_Y","PRESSURE"],
                "input_file_name"      : "reference_couette_flow_test",
                "model_part_name"      : "FluidModelPart.FluidParts_Fluid",
                "tolerance"            : 0.0,
                "relative_tolerance"   : 0.0,
                "time_frequency"       : 99.0
            }
        }""")
        input_file_name = json_check_settings["Parameters"]["input_file_name"].GetString()
        input_file_name += "_" + self.element_type + ".json"
        json_check_settings["Parameters"]["input_file_name"].SetString(input_file_name)
        json_check_settings["Parameters"]["tolerance"].SetDouble(self.check_tolerance)
        json_check_settings["Parameters"]["relative_tolerance"].SetDouble(self.check_relative_tolerance)
        self.parameters["processes"]["json_check_process_list"].Append(json_check_settings)

if __name__ == '__main__':
    KratosUnittest.main()
