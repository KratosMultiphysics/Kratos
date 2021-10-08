import KratosMultiphysics
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.FluidDynamicsApplication import initialize_with_compressible_potential_solution_process
from KratosMultiphysics import KratosUnittest
class InitializeWithCompressiblePotentialSolutionProcessTest(KratosUnittest.TestCase):
    def setUp(self):
        pass

    @KratosUnittest.skipIfApplicationsNotAvailable("CompressiblePotentialFlowApplication")
    def testAnalysisParameters(self):
        work_folder = "Square"
        work_file = "square"
        with KratosUnittest.WorkFolderScope(work_folder, __file__):
            model = KratosMultiphysics.Model()
            mpart = model.CreateModelPart("main_model_part")
            mpart.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
            ReadModelPart(work_file, mpart)

            process = initialize_with_compressible_potential_solution_process.Factory(self._GetSettings(), model)

            obtained = process._GenerateAnalysisparameters(process.settings, 0.0)
            expected = self._GetExpectedAnalysisParameters()

            obtained.RecursivelyValidateDefaults(expected)

            self.assertTrue(
                obtained.IsEquivalentTo(expected),
                "Parameters do not match. Obtained:\n\n{}\n\nExpected:\n\n{}".
                    format(obtained.PrettyPrintJsonString(), expected.PrettyPrintJsonString())
            )

    @KratosUnittest.skipIfApplicationsNotAvailable("CompressiblePotentialFlowApplication")
    def testProcess(self):
        "Ensures the whole process functions when given expected parameters"
        work_folder = "Square"
        work_file = "square"
        with KratosUnittest.WorkFolderScope(work_folder, __file__):
            model = KratosMultiphysics.Model()
            mpart = model.CreateModelPart("main_model_part")

            mpart.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
            mpart.SetBufferSize(2)
            mpart.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
            mpart.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENTUM)
            mpart.AddNodalSolutionStepVariable(KratosMultiphysics.TOTAL_ENERGY)
            mpart.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
            mpart.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
            ReadModelPart(work_file, mpart)

            process = initialize_with_compressible_potential_solution_process.Factory(self._GetSettings(), model)
            process.ExecuteInitialize()
            process.ExecuteBeforeSolutionLoop()

            expected_v = 0.8 * 344.31
            expected_rho = 1.19659
            expected_e = 298706

            for node in mpart.Nodes:
                self.assertAlmostEqual(
                    node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X),
                    expected_v,
                    places=7)
                self.assertAlmostEqual(
                    node.GetSolutionStepValue(KratosMultiphysics.DENSITY),
                    expected_rho,
                    places=5)
                self.assertAlmostEqual(
                    node.GetSolutionStepValue(KratosMultiphysics.MOMENTUM_X),
                    expected_v * expected_rho,
                    places=2)
                self.assertAlmostEqual(
                    node.GetSolutionStepValue(KratosMultiphysics.TOTAL_ENERGY),
                    expected_e,
                    places=0) # Large value -> Lower absolute accuracy


    @classmethod
    def _GetSettings(cls):
        """Returns settings for M=0.8 @ 1bar and 293.15K"""
        return KratosMultiphysics.Parameters("""
            {
                "Parameters" : {
                    "model_part_name": "main_model_part",
                    "volume_model_part_name": "Fluid",
                    "skin_parts": ["Inlet","Outlet"],
                    "properties" : {
                        "c_v" : 722.14,
                        "gamma" : 1.4,
                        "free_stream_density" : 1.19659,
                        "free_stream_momentum" : 329.598,
                        "free_stream_energy" : 298706.0
                    },
                    "boundary_conditions_process_list": [
                        {
                            "python_module" : "apply_far_field_process",
                            "kratos_module" : "KratosMultiphysics.CompressiblePotentialFlowApplication",
                            "process_name"  : "FarFieldProcess",
                            "Parameters"    : {
                                "model_part_name" : "potential_analysis_model_part.Inlet",
                                "angle_of_attack" : 10.0,
                                "mach_infinity"   : 0.8,
                                "speed_of_sound"  : 344.31
                            }
                        },{
                            "python_module" : "apply_far_field_process",
                            "kratos_module" : "KratosMultiphysics.CompressiblePotentialFlowApplication",
                            "process_name"  : "FarFieldProcess",
                            "Parameters"    : {
                                "model_part_name" : "potential_analysis_model_part.Outlet",
                                "angle_of_attack" : 0.0,
                                "mach_infinity"   : 0.8,
                                "speed_of_sound"  : 344.31
                            }
                        }
                    ]
                }
            }
        """)

    def _GetExpectedAnalysisParameters(self):
        return KratosMultiphysics.Parameters("""
            {
                "output_processes": {},
                "problem_data": {
                    "echo_level": 0,
                    "end_time": 0.0,
                    "parallel_type": "OpenMP",
                    "problem_name": "potential_analysis_model_part",
                    "start_time": -1.0
                },
                "processes": {
                    "auxiliar_process_list": [],
                    "boundary_conditions_process_list": [
                        {
                            "Parameters": {
                                "angle_of_attack": 10.0,
                                "mach_infinity": 0.8,
                                "model_part_name": "potential_analysis_model_part.Inlet",
                                "speed_of_sound": 344.31
                            },
                            "kratos_module": "KratosMultiphysics.CompressiblePotentialFlowApplication",
                            "process_name": "FarFieldProcess",
                            "python_module": "apply_far_field_process"
                        },
                        {
                            "Parameters": {
                                "angle_of_attack": 0.0,
                                "mach_infinity": 0.8,
                                "model_part_name": "potential_analysis_model_part.Outlet",
                                "speed_of_sound": 344.31
                            },
                            "kratos_module": "KratosMultiphysics.CompressiblePotentialFlowApplication",
                            "process_name": "FarFieldProcess",
                            "python_module": "apply_far_field_process"
                        }
                    ]
                },
                "solver_settings": {
                    "domain_size": 2,
                    "echo_level": 0,
                    "maximum_iterations": 10,
                    "model_import_settings": {
                        "input_type": "use_input_model_part"
                    },
                    "formulation": {
                        "element_type" : "compressible"
                    },
                    "model_part_name": "potential_analysis_model_part",
                    "no_skin_parts": [],
                    "reform_dofs_at_each_step": false,
                    "skin_parts": [
                        "Inlet",
                        "Outlet"
                    ],
                    "solver_type": "potential_flow",
                    "volume_model_part_name": "Fluid"
                }
            }

        """)

if __name__ == '__main__':
    KratosUnittest.main()