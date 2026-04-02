# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.LinearSolversApplication as KratosLinearSolvers

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver


def CreateSolver(model, custom_settings):
    return HarmonicAnalysisSolver(model, custom_settings)


class HarmonicAnalysisSolver(MechanicalSolver):
    """The structural mechanics harmonic analysis solver.

    Supports:
    - modal harmonic analysis
    - direct harmonic analysis
    """

    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)

        harmonic_type = self.settings["harmonic_analysis_type"].GetString()
        if harmonic_type not in ["modal", "direct"]:
            err_msg = f'Unknown "harmonic_analysis_type": "{harmonic_type}"\n'
            err_msg += 'Available options are: "modal", "direct"'
            raise Exception(err_msg)

        self._ValidateCustomSettings()

        KratosMultiphysics.Logger.PrintInfo(
            "::[HarmonicAnalysisSolver]:: ",
            f'Construction finished, harmonic_analysis_type = "{harmonic_type}"'
        )

    def _ValidateCustomSettings(self):
        self.settings["modal_harmonic_analysis_settings"].ValidateAndAssignDefaults(
            KratosMultiphysics.Parameters(r'''{
                "use_effective_material_damping" : false
            }''')
        )

        self.settings["direct_harmonic_analysis_settings"].ValidateAndAssignDefaults(
            KratosMultiphysics.Parameters(r'''{
                "mass_matrix_diagonal_value"      : 1.0,
                "stiffness_matrix_diagonal_value" : 1.0,
                "damping_matrix_diagonal_value"   : 1.0,
                "reform_dof_set_at_each_step"     : false,
                "assemble_damping_matrix"         : false,
                "real_load_sub_model_part"        : "",
                "imaginary_load_sub_model_part" : "",
                "complex_linear_solver_settings"  : {
                    "solver_type" : "pardiso_lu_complex"
                }
            }''')
        )
    
    def AddVariables(self):
        # First add the standard structural mechanics variables
        super().AddVariables()

        # Then add harmonic-analysis-specific variables
        self.main_model_part.AddNodalSolutionStepVariable(
            StructuralMechanicsApplication.DISPLACEMENT_IMAGINARY
        )
        self.main_model_part.AddNodalSolutionStepVariable(
            StructuralMechanicsApplication.POINT_LOAD_IMAGINARY
        )

        KratosMultiphysics.Logger.PrintInfo(
            "::[HarmonicMechanicalSolver]:: ",
            "Harmonic variables ADDED"
        )

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters(r'''{
            "scheme_type" : "dynamic",
            "harmonic_analysis_type" : "modal",
            "modal_harmonic_analysis_settings" : {
                "use_effective_material_damping" : false
            },
            "direct_harmonic_analysis_settings" : {
                "mass_matrix_diagonal_value"      : 1.0,
                "stiffness_matrix_diagonal_value" : 1.0,
                "damping_matrix_diagonal_value"   : 1.0,
                "reform_dof_set_at_each_step"     : false,
                "assemble_damping_matrix"         : false,
                "real_load_sub_model_part"        : "",
                "imaginary_load_sub_model_part" : "",
                "complex_linear_solver_settings"  : {
                    "solver_type" : "pardiso_lu_complex"
                }
            }
        }''')
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    #### Private functions ####

    def _CreateScheme(self):
        """Create the scheme to construct the global operators."""
        scheme_type = self.settings["scheme_type"].GetString()

        if scheme_type == "dynamic":
            return StructuralMechanicsApplication.EigensolverDynamicScheme()
        else:
            err_msg = f'The requested scheme type "{scheme_type}" is not available!\n'
            err_msg += 'Available options are: "dynamic"'
            raise Exception(err_msg)

    def _CreateLinearSolver(self):
        """Create the real linear solver used by the standard builder-and-solver.

        Important:
        - The direct harmonic strategy uses its own COMPLEX linear solver.
        - The standard StructuralMechanics builder-and-solver still expects a
          regular real Kratos.LinearSolver-compatible object.
        """
        harmonic_type = self.settings["harmonic_analysis_type"].GetString()

        if harmonic_type == "modal":
            # Same behavior as before: modal strategy does not need a standard solve.
            return KratosMultiphysics.LinearSolver()

        elif harmonic_type == "direct":
            # Keep the builder-and-solver infrastructure real-valued.
            # The complex solver is created separately in _CreateComplexLinearSolver().
            return KratosMultiphysics.LinearSolver()

        else:
            err_msg = f'Unknown "harmonic_analysis_type": "{harmonic_type}"'
            raise Exception(err_msg)

    def _CreateComplexLinearSolver(self):
        """Create the complex linear solver used ONLY by the direct harmonic strategy.
        """
        KratosMultiphysics.Logger.PrintInfo(
            "::[HarmonicAnalysisSolver]::",
            self.settings.PrettyPrintJsonString()
        )
        complex_solver_settings = self.settings["direct_harmonic_analysis_settings"]["complex_linear_solver_settings"]

        factory = KratosMultiphysics.ComplexLinearSolverFactory()
        return factory.Create(complex_solver_settings)

    def _CreateSolutionStrategy(self):
        scheme = self._GetScheme()
        builder_and_solver = self._GetBuilderAndSolver()
        computing_model_part = self.GetComputingModelPart()

        harmonic_type = self.settings["harmonic_analysis_type"].GetString()

        if harmonic_type == "modal":
            return StructuralMechanicsApplication.ModalHarmonicAnalysisStrategy(
                computing_model_part,
                scheme,
                builder_and_solver,
                self.settings["modal_harmonic_analysis_settings"]["use_effective_material_damping"].GetBool()
            )

        elif harmonic_type == "direct":
            direct_settings = self.settings["direct_harmonic_analysis_settings"]
            complex_linear_solver = self._CreateComplexLinearSolver()

            return StructuralMechanicsApplication.DirectHarmonicAnalysisStrategy(
                computing_model_part,
                scheme,
                builder_and_solver,
                complex_linear_solver,
                direct_settings["mass_matrix_diagonal_value"].GetDouble(),
                direct_settings["stiffness_matrix_diagonal_value"].GetDouble(),
                direct_settings["damping_matrix_diagonal_value"].GetDouble(),
                direct_settings["reform_dof_set_at_each_step"].GetBool(),
                direct_settings["assemble_damping_matrix"].GetBool(),
                direct_settings["real_load_sub_model_part"].GetString(),
                direct_settings["imaginary_load_sub_model_part"].GetString()
            )

        else:
            err_msg = f'Unknown "harmonic_analysis_type": "{harmonic_type}"'
            raise Exception(err_msg)