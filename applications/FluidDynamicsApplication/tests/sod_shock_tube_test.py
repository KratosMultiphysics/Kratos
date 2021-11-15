# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

class SodShockTubeTest(KratosUnittest.TestCase):
    def testSodShockTubeExplicitASGS(self):
        self.solver_type = "CompressibleExplicit"
        self.use_oss = False
        self.shock_capturing = False
        self._CustomizeTestSettings()
        self._RunSodShockTubeTest()

    def testSodShockTubeExplicitASGSShockCapturing(self):
        self.solver_type = "CompressibleExplicit"
        self.use_oss = False
        self.shock_capturing = True
        self._CustomizeTestSettings()
        self._RunSodShockTubeTest()

    def testSodShockTubeExplicitOSS(self):
        self.solver_type = "CompressibleExplicit"
        self.use_oss = True
        self.shock_capturing = False
        self._CustomizeTestSettings()
        self._RunSodShockTubeTest()

    def testSodShockTubeExplicitOSSShockCapturing(self):
        self.solver_type = "CompressibleExplicit"
        self.use_oss = True
        self.shock_capturing = True
        self._CustomizeTestSettings()
        self._RunSodShockTubeTest()

    def setUp(self):
        self.print_output = False
        self.print_reference_values = False
        self.check_absolute_tolerance = 1.0e-8
        self.check_relative_tolerance = 1.0e-10
        self.work_folder = "sod_shock_tube_test"
        settings_filename = "ProjectParameters.json"

        # Read the simulation settings
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            with open(settings_filename,'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            KratosUtilities.DeleteFileIfExisting('sod_shock_tube_geom_coarse.time')

    def _RunSodShockTubeTest(self):
        # Create the test simulation
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            self.model = KratosMultiphysics.Model()
            simulation = FluidDynamicsAnalysis(self.model, self.parameters)
            simulation.Run()

    def _CustomizeTestSettings(self):
        # Customize simulation settings
        self.parameters["solver_settings"]["solver_type"].SetString(self.solver_type)
        self.parameters["solver_settings"]["use_oss"].SetBool(self.use_oss)
        self.parameters["solver_settings"]["shock_capturing"].SetBool(self.shock_capturing)

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
                "model_part_name"        : "FluidModelPart",
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
                        "nodal_results"               : ["DENSITY","MOMENTUM","TOTAL_ENERGY"],
                        "gauss_point_results"         : ["SHOCK_SENSOR","THERMAL_SENSOR","SHEAR_SENSOR"],
                        "nodal_nonhistorical_results" : ["ARTIFICIAL_BULK_VISCOSITY","ARTIFICIAL_CONDUCTIVITY","ARTIFICIAL_DYNAMIC_VISCOSITY"]
                    },
                    "point_data_configuration"  : []
                }
            }
        }""")
        output_name = "sod_shock_tube{0}{1}{2}".format(
            "_explicit" if self.solver_type == "CompressibleExplicit" else "_implicit",
            "_ASGS" if self.use_oss == False else "_OSS",
            "_SC" if self.shock_capturing else "")
        gid_output_settings["Parameters"]["output_name"].SetString(output_name)
        self.parameters["output_processes"]["gid_output"].Append(gid_output_settings)

    def _AddReferenceValuesOutput(self):
        json_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "json_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "JsonOutputProcess",
            "Parameters"    : {
                "output_variables" : ["DENSITY","MOMENTUM_X","MOMENTUM_Y","TOTAL_ENERGY"],
                "output_file_name" : "TO_BE_DEFINED",
                "model_part_name"  : "FluidModelPart.FluidParts_Fluid",
                "time_frequency"   : 0.025
            }
        }""")
        output_file_name = "sod_shock_tube{0}{1}{2}_results.json".format(
            "_explicit" if self.solver_type == "CompressibleExplicit" else "_implicit",
            "_ASGS" if self.use_oss == False else "_OSS",
            "_SC" if self.shock_capturing else "")
        json_output_settings["Parameters"]["output_file_name"].SetString(output_file_name)
        self.parameters["processes"]["json_check_process_list"].Append(json_output_settings)

    def _AddReferenceValuesCheck(self):
        json_check_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "from_json_check_result_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "FromJsonCheckResultProcess",
            "Parameters"    : {
                "check_variables"      : ["DENSITY","MOMENTUM_X","MOMENTUM_Y","TOTAL_ENERGY"],
                "input_file_name"      : "TO_BE_DEFINED",
                "model_part_name"      : "FluidModelPart.FluidParts_Fluid",
                "tolerance"            : 0.0,
                "relative_tolerance"   : 0.0,
                "time_frequency"       : 0.025
            }
        }""")
        input_file_name = "sod_shock_tube{0}{1}{2}_results.json".format(
            "_explicit" if self.solver_type == "CompressibleExplicit" else "_implicit",
            "_ASGS" if self.use_oss == False else "_OSS",
            "_SC" if self.shock_capturing else "")
        json_check_settings["Parameters"]["input_file_name"].SetString(input_file_name)
        json_check_settings["Parameters"]["tolerance"].SetDouble(self.check_absolute_tolerance)
        json_check_settings["Parameters"]["relative_tolerance"].SetDouble(self.check_relative_tolerance)
        self.parameters["processes"]["json_check_process_list"].Append(json_check_settings)

if __name__ == '__main__':
    KratosUnittest.main()
    # test = SodShockTubeTest()
    # test.testSodShockTubeExplicitASGS()
    # test.testSodShockTubeExplicitASGSShockCapturing()
    # test.testSodShockTubeExplicitOSS()
    # test.testSodShockTubeExplicitOSSShockCapturing()