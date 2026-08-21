# --- Core Imports ---
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.analysis_stage import AnalysisStage

# --- Structural Mechanics Imports ---
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

# --- STD Imports ---
import pathlib


class PMGBuilderAndSolverSingleGridTest(UnitTest.TestCase):

    @staticmethod
    def __ApplySetting(parameters: KratosMultiphysics.Parameters,
                       key: str,
                       value: KratosMultiphysics.Parameters) -> None:
        if parameters.Has(key):
            parameters.RemoveValue(key)
        parameters.AddValue(key, value)

    def Run(self, directory: pathlib.Path) -> None:
        # Define constraint assemblers to test.
        constraint_assemblers: dict[str,KratosMultiphysics.Parameters] = {
            "master_slave" : KratosMultiphysics.Parameters(
                R"""{
                    "method" : "master_slave"
                }"""),
            "augmented_lagrange" : KratosMultiphysics.Parameters(
                R"""{
                    "method" : "augmented_lagrange",
                    "max_iterations" : 5e1,
                    "tolerance" : 1e-12,
                    "penalty_factor" : "1e3 * max"
                }""")
        }

        with UnitTest.WorkFolderScope(".", __file__):
            config_files: list[pathlib.Path]
            if directory.is_file():
                config_files = [directory]
            else:
                config_files = [p for p in directory.rglob("*_parameters.json")]
            if not config_files:
                raise RuntimeError(f"no tests found in {directory}")
            for parameters_file in config_files:
                for constraint_imposition_name, constraint_imposition_settings in constraint_assemblers.items():
                    with self.subTest(f"imposition: {constraint_imposition_name} configuration: {parameters_file}"):

                        parameters: KratosMultiphysics.Parameters
                        with open(parameters_file, "r") as file:
                            parameters = KratosMultiphysics.Parameters(file.read())

                        # Assemble the input configuration for PMGBuilderAndSolver
                        pmg_parameters: KratosMultiphysics.Parameters = KratosMultiphysics.Parameters(
                            R"""{
                                "type" : "p_multigrid",
                                "advanced_settings" : {
                                    "max_iterations" : 1,
                                    "verbosity" : 0,
                                    "linear_solver_settings" : {
                                        "solver_type" : "skyline_lu_factorization"
                                    }
                                }
                            }""")
                        pmg_parameters["advanced_settings"].AddValue(
                            "constraint_imposition_settings",
                            constraint_imposition_settings)

                        # Set PMGBuilderAndSolver in the input configuration.
                        if parameters["solver_settings"].Has("block_builder"):
                            parameters["solver_settings"].RemoveValue("block_builder")
                        self.__ApplySetting(
                            parameters["solver_settings"],
                            "builder_and_solver_settings",
                            pmg_parameters)

                        # Reduce analysis verbosity.
                        self.__ApplySetting(
                            parameters["problem_data"],
                            "echo_level",
                            KratosMultiphysics.Parameters("0"))
                        self.__ApplySetting(
                            parameters["solver_settings"],
                            "echo_level",
                            KratosMultiphysics.Parameters("0"))

                        model: KratosMultiphysics.Model = KratosMultiphysics.Model()
                        analysis: AnalysisStage = StructuralMechanicsAnalysis(model, parameters)
                        analysis.Run()

    def test_Execution(self) -> None:
        for directory in [
            pathlib.Path("reissner_mindlin_shells"),
            pathlib.Path("mixed_u_E_test"),
            pathlib.Path("LinearTruss2D"),
            pathlib.Path("LinearTruss3D"),
            pathlib.Path("TLTruss3D"),
            pathlib.Path("TimoshenkoBeams"),
            pathlib.Path("automated_initial_variable_process_test"),
            #pathlib.Path("patch_test"), # <== fucked, shifted boundary tests fail with standard settings too
            #pathlib.Path("sprism_test"),
            pathlib.Path("explicit_tests"),
            pathlib.Path("formfinding_test"),
            pathlib.Path("membrane_test"),
            pathlib.Path("truss_test") / "dynamic_3D2NTruss_test_parameters.json",
            pathlib.Path("truss_test") / "linear_3D2NTruss_plastic_compression_test_parameters.json",
            pathlib.Path("truss_test") / "linear_3D2NTruss_plastic_tension_test_parameters.json",
            pathlib.Path("truss_test") / "linear_3D2NTruss_test_parameters.json",
            pathlib.Path("truss_test") / "nonlinear_3D2NTruss_plastic_snapthrough_test_parameters.json",
            pathlib.Path("truss_test") / "nonlinear_3D2NTruss_plastic_tension_test_parameters.json",
            pathlib.Path("truss_test") / "nonlinear_3D2NTruss_test_parameters.json",
            pathlib.Path("beam_test"),
            #pathlib.Path("InitialStateElasticity"), # <== fucked, InitialStateElasticity/initial_strain_cr_beam3D2N_test_parameters.json points to an invalid mesh path
            pathlib.Path("shell_test") / "Shell_T3_isotropic_linear_static_struct_scordelis_lo_roof_parameters.json",
            pathlib.Path("shell_test") / "Shell_T3andQ4_linear_static_struct_scordelis_lo_roof_parameters.json",
            pathlib.Path("shell_test") / "Shell_T3andQ4_linear_static_struct_pinched_cylinder_parameters.json",
            pathlib.Path("shell_test") / "Shell_T3andQ4_linear_static_struct_pinched_hemisphere_parameters.json",
            pathlib.Path("shell_test") / "Shell_T3andQ4_linear_static_struct_clamped_cylinder_orthotropic_parameters.json",
            pathlib.Path("shell_test") / "Shell_T3andQ4_nonlinear_static_struct_hinged_cyl_roof_snapthrough_parameters.json",
            pathlib.Path("shell_test") / "Shell_T3andQ4_nonlinear_static_unstruct_hinged_cyl_roof_snapthrough_parameters.json",
            #pathlib.Path("shell_test") / "Shell_T3andQ4_nonlinear_dynamic_struct_oscillating_plate_parameters.json", # <== broken results path
            #pathlib.Path("shell_test") / "Shell_T3andQ4_nonlinear_dynamic_struct_oscillating_plate_lumped_parameters.json", # <== Value type for "USE_LUMPED_MASS_MATRIX" not defined
            pathlib.Path("rigid_test"),
            pathlib.Path("solid_2p5d_test"),
            ]:
                self.Run(directory)


if __name__ == "__main__":
    UnitTest.main()
