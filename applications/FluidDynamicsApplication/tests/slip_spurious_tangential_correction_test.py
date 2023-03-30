# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

class SlipSpuriousTangentialCorrectionTest(KratosUnittest.TestCase):
    def setUp(self):
        self.print_output = False
        self.check_tolerance = 1.0e-6
        self.check_relative_tolerance = 1.0e-8
        self.print_reference_values = False
        self.work_folder = "slip_spurious_tangential_correction_test"
        settings_filename = "ProjectParameters.json"

        # Read the simulation settings
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            with open(settings_filename,'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())

        # Customize the test settings
        self.__CustomizeTestSettings()

    def testSpuriousTangentialCorrectionTest2D(self):
        # Create the test simulation
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            self.model = KratosMultiphysics.Model()
            simulation = FluidDynamicsAnalysis(self.model, self.parameters)
            simulation.Run()

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            KratosUtilities.DeleteFileIfExisting('slip_spurious_tangential_correction.time')

    def __CustomizeTestSettings(self):
        # If required, add the output process to the test settings
        if self.print_output:
            self.__AddOutput()

        # If required, add the reference values output process to the test settings
        if self.print_reference_values:
            self.__AddReferenceValuesOutput()
        else:
            self.__AddReferenceValuesCheck()

    def __AddOutput(self):
        gid_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "help"          : "This process writes postprocessing files for GiD",
            "Parameters"    : {
                "model_part_name"        : "FluidModelPart.fluid_computational_model_part",
                "output_name"            : "slip_spurious_tangential_correction",
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
        self.parameters["output_processes"]["gid_output"].Append(gid_output_settings)

    def __AddReferenceValuesOutput(self):
        json_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "json_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "JsonOutputProcess",
            "Parameters"    : {
                "output_variables" : ["VELOCITY","PRESSURE"],
                "output_file_name" : "reference_slip_spurious_tangential_correction.json",
                "model_part_name"  : "FluidModelPart.FluidParts_Fluid",
                "time_frequency"   : 0.99
            }
        }""")
        self.parameters["processes"]["json_check_process_list"].Append(json_output_settings)

    def __AddReferenceValuesCheck(self):
        json_check_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "from_json_check_result_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "FromJsonCheckResultProcess",
            "Parameters"    : {
                "check_variables"      : ["VELOCITY_X","VELOCITY_Y","PRESSURE"],
                "input_file_name"      : "reference_slip_spurious_tangential_correction.json",
                "model_part_name"      : "FluidModelPart.FluidParts_Fluid",
                "tolerance"            : 0.0,
                "relative_tolerance"   : 0.0,
                "time_frequency"       : 99.0
            }
        }""")
        json_check_settings["Parameters"]["tolerance"].SetDouble(self.check_tolerance)
        json_check_settings["Parameters"]["relative_tolerance"].SetDouble(self.check_relative_tolerance)
        self.parameters["processes"]["json_check_process_list"].Append(json_check_settings)

if __name__ == '__main__':
    KratosUnittest.main()
