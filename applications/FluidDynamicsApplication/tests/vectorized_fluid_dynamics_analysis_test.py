# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities
import KratosMultiphysics.FluidDynamicsApplication.cfd_utils as cfd_utils
from KratosMultiphysics.FluidDynamicsApplication.vectorized_cfd_stage import VectorizedCFDStage

# Check if cupy is available for GPU tests
try:
    import cupy
    HAS_CUPY = True
except ImportError:
    HAS_CUPY = False
    print("Skipping tests that require CuPy.")

class _VectorizedFluidDynamicsAnalysisTest(KratosUnittest.TestCase):

    def setUp(self):
        self.print_output = False
        self.print_reference_values = False
        self.check_absolute_tolerance = 1.0e-4
        self.check_relative_tolerance = 1.0e-6
        self.work_folder = "vectorized_fluid_dynamics_analysis_test"
        settings_filename = "ProjectParameters.json"

        # Read the simulation settings
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            with open(settings_filename,'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())

    def tearDown(self):
        # Reset the backend between tests
        cfd_utils.reset()

        # Wipe the auxiliary mdpa reading files
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            KratosUtilities.DeleteFileIfExisting('lid_driven_cavity_10x10.time')

    def _RunVectorizedAnalysis(self):
        # Create the test simulation
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            self.model = KratosMultiphysics.Model()
            simulation = VectorizedCFDStage(self.model, self.parameters)
            simulation.Run()

    def _CustomizeTestSettings(self):
        # Customize simulation settings
        self.parameters["problem_data"]["parallel_type"].SetString(self.parallel_type)
        self.parameters["problem_data"]["precision"].SetString(self.precision)

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
            "Parameters"    : {
                "model_part_name"        : "FluidModelPart",
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
                        "output_interval"             : 1.0,
                        "body_output"                 : true,
                        "node_output"                 : false,
                        "skin_output"                 : false,
                        "plane_output"                : [],
                        "nodal_results"               : ["VELOCITY","PRESSURE"],
                        "nodal_nonhistorical_results" : [],
                        "gauss_point_results"         : []
                    },
                    "point_data_configuration"  : []
                },
                "output_name"            : ""
            }
        }""")
        output_name = f"lid_driven_cavity_10x10_{self.parallel_type}_{self.precision}"
        gid_output_settings["Parameters"]["output_name"].SetString(output_name)
        self.parameters["output_processes"]["gid_output"].Append(gid_output_settings)

    def _AddReferenceValuesOutput(self):
        json_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "json_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "JsonOutputProcess",
            "Parameters"    : {
                "output_variables" : ["PRESSURE","VELOCITY_X","VELOCITY_Y"],
                "output_file_name" : "TO_BE_DEFINED",
                "model_part_name"  : "FluidModelPart",
                "time_frequency"   : 0.025
            }
        }""")
        output_file_name = f"lid_driven_cavity_10x10_{self.parallel_type}_{self.precision}_results.json"
        json_output_settings["Parameters"]["output_file_name"].SetString(output_file_name)
        self.parameters["processes"]["json_check_process_list"].Append(json_output_settings)

    def _AddReferenceValuesCheck(self):
        json_check_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "from_json_check_result_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "FromJsonCheckResultProcess",
            "Parameters"    : {
                "check_variables"      : ["PRESSURE","VELOCITY_X","VELOCITY_Y"],
                "input_file_name"      : "TO_BE_DEFINED",
                "model_part_name"      : "FluidModelPart",
                "tolerance"            : 0.0,
                "relative_tolerance"   : 0.0,
                "time_frequency"       : 0.025
            }
        }""")
        input_file_name = f"lid_driven_cavity_10x10_{self.parallel_type}_{self.precision}_results.json"
        json_check_settings["Parameters"]["input_file_name"].SetString(input_file_name)
        json_check_settings["Parameters"]["tolerance"].SetDouble(self.check_absolute_tolerance)
        json_check_settings["Parameters"]["relative_tolerance"].SetDouble(self.check_relative_tolerance)
        self.parameters["processes"]["json_check_process_list"].Append(json_check_settings)

@KratosUnittest.skipUnless(HAS_CUPY, "CuPy not available")
class VectorizedFluidDynamicsAnalysisTestGPU(_VectorizedFluidDynamicsAnalysisTest):
    def test_gpu_float32(self):
        self.parallel_type = "gpu"
        self.precision = "float32"
        self._CustomizeTestSettings()
        self._RunVectorizedAnalysis()

@KratosUnittest.skipUnless(HAS_CUPY, "CuPy not available")
class VectorizedFluidDynamicsAnalysisTestOpenMP(_VectorizedFluidDynamicsAnalysisTest):
    def test_open_mp_float32(self):
        self.parallel_type = "open_mp"
        self.precision = "float32"
        self._CustomizeTestSettings()
        self._RunVectorizedAnalysis()

if __name__ == '__main__':
    KratosUnittest.main()
    # test = SodShockTubeTest()
    # test.testSodShockTubeExplicitASGS()
    # test.testSodShockTubeExplicitASGSShockCapturing()
    # test.testSodShockTubeExplicitOSS()
    # test.testSodShockTubeExplicitOSSShockCapturing()