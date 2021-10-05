import KratosMultiphysics
import KratosMultiphysics.KratosUnittest
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.FluidDynamicsApplication import initialize_with_compressible_potential_solution_process

from KratosMultiphysics import KratosUnittest

class InitializeWithCompressiblePotentialSolutionProcessTest(KratosUnittest.TestCase):
    def setUp(self):
        pass

    # @KratosUnittest.skipIfApplicationsNotAvailable("CompressiblePotentialFlowApplication")
    # def test_AnalysisParameters(self):
    #     work_folder = "Square"
    #     work_file = "square"
    #     with KratosUnittest.WorkFolderScope(work_folder, __file__):
    #         model = KratosMultiphysics.Model()
    #         mpart = model.CreateModelPart("main_model_part")
    #         mpart.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
    #         ReadModelPart(work_file, mpart)

    #         process = initialize_with_compressible_potential_solution_process.Factory(self._GetSettings(), model)

    #         obtained = process.analysis_parameters
    #         expected = self._GetExpectedAnalysisParameters()

    #         obtained.RecursivelyValidateDefaults(expected)

    #         self.assertTrue(
    #             obtained.IsEquivalentTo(expected),
    #             "Parameters do not match. Obtained:\n\n{}\n\nExpected:\n\n{}".
    #                 format(obtained.PrettyPrintJsonString(), expected.PrettyPrintJsonString())
    #         )

    @KratosUnittest.skipIfApplicationsNotAvailable("CompressiblePotentialFlowApplication")
    def test_Process(self):
        work_folder = "Square"
        work_file = "square"
        with KratosUnittest.WorkFolderScope(work_folder, __file__):
            model = KratosMultiphysics.Model()
            mpart = model.CreateModelPart("main_model_part")
            mpart.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
            ReadModelPart(work_file, mpart)

            process = initialize_with_compressible_potential_solution_process.Factory(self._GetSettings(), model)
            process.ExecuteInitialize()



    @classmethod
    def _GetSettings(cls):
        """Returns settings for M=0.8 @ 1bar and 293.15K"""
        return KratosMultiphysics.Parameters("""
            {
                "Parameters" : {
                    "model_part_name": "main_model_part",
                    "volume_model_part_name": "model_part_potential_analysis.Fluid",
                    "skin_parts": ["Inlet","Outlet","Walls"],
                    "materials_filename": "fluid_materials.json",
                    "free_stream_density": 1.19659,
                    "free_stream_momentum": 329.598,
                    "free_stream_energy": 298706.0,
                    "boundary_conditions_process_list": [
                        {
                            "python_module": "apply_far_field_process",
                            "kratos_module": "KratosMultiphysics.CompressiblePotentialFlowApplication",
                            "process_name": "FarFieldProcess",
                            "Parameters":
                            {
                                "model_part_name": "model_part_potential_analysis.Inlet",
                                "angle_of_attack": 0.0,
                                "mach_infinity": 0.8,
                                "speed_of_sound": 344.31
                            }
                        },
                        {
                            "python_module": "apply_far_field_process",
                            "kratos_module": "KratosMultiphysics.CompressiblePotentialFlowApplication",
                            "process_name": "FarFieldProcess",
                            "Parameters":
                            {
                                "model_part_name": "model_part_potential_analysis.Outlet",
                                "angle_of_attack": 0.0,
                                "mach_infinity": 0.8,
                                "speed_of_sound": 344.31
                            }
                        },
                        {
                            "python_module": "define_wake_process_2d",
                            "kratos_module": "KratosMultiphysics.CompressiblePotentialFlowApplication",
                            "process_name": "DefineWakeProcess2D",
                            "Parameters":
                            {
                                "model_part_name": "model_part_potential_analysis.Walls"
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
                "end_time": 1,
                "parallel_type": "OpenMP",
                "problem_name": "internal_potential_solver",
                "start_time": 0.0
            },
            "processes": {
                "auxiliar_process_list": [],
                "boundary_conditions_process_list": [
                    {
                        "Parameters": {
                            "angle_of_attack": 0.0,
                            "mach_infinity": 0.8,
                            "model_part_name": "FluidModelPart.AutomaticInlet2D_Left",
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
                            "model_part_name": "FluidModelPart.Outlet2D_Right",
                            "speed_of_sound": 344.31
                        },
                        "kratos_module": "KratosMultiphysics.CompressiblePotentialFlowApplication",
                        "process_name": "FarFieldProcess",
                        "python_module": "apply_far_field_process"
                    },
                    {
                        "Parameters": {
                            "model_part_name": "FluidModelPart.NoSlip2D_Aerofoil"
                        },
                        "kratos_module": "KratosMultiphysics.CompressiblePotentialFlowApplication",
                        "process_name": "DefineWakeProcess2D",
                        "python_module": "define_wake_process_2d"
                    }
                ]
            },
            "solver_settings": {
                "domain_size": 2,
                "echo_level": 0,
                "formulation": {
                    "element_type": "compressible"
                },
                "material_import_settings": {
                    "materials_filename": "fluid_materials.json"
                },
                "maximum_iterations": 10,
                "model_import_settings": {
                    "input_type": "use_input_model_part"
                },
                "model_part_name": "main_model_part",
                "no_skin_parts": [],
                "reform_dofs_at_each_step": false,
                "skin_parts": [
                    "Inlet",
                    "Outlet",
                    "Walls"
                ],
                "solver_type": "potential_flow",
                "volume_model_part_name": "Fluid"
            }
        }
        """)

if __name__ == '__main__':
    KratosUnittest.main()