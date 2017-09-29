from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import structural_mechanics_solver

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(main_model_part, custom_settings):
    return HarmonicAnalysisSolver(main_model_part, custom_settings)


class HarmonicAnalysisSolver(structural_mechanics_solver.MechanicalSolver):
    """The structural mechanics harmonic analysis solver.

    This class creates the mechanical solvers for the harmonic analysis.
    It currently supports the Feast solver.

    Public member variables:
    harmonic_analysis_settings -- settings for the eigenvalue solvers.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        # Set defaults and validate custom settings.
        harmonic_analysis_settings = KratosMultiphysics.Parameters("""
        {
            "harmonic_analysis_settings" : {
                "system_damping" : 0.0,
                "rayleigh_alpha" : 0.0,
                "rayleigh_beta"  : 0.0
            }
        }
        """)
        self.validate_and_transfer_matching_settings(custom_settings, harmonic_analysis_settings)
        self.harmonic_analysis_settings = harmonic_analysis_settings["harmonic_analysis_settings"]
        # Validate the remaining settings in the base class.
        if not custom_settings.Has("scheme_type"): # Override defaults in the base class.
            custom_settings.AddEmptyValue("scheme_type")
            custom_settings["scheme_type"].SetString("Dynamic")
        
        # Construct the base solver.
        super(HarmonicAnalysisSolver, self).__init__(main_model_part, custom_settings)
        print("::[HarmonicAnalysisSolver]:: Construction finished")

    #### Private functions ####

    def _create_solution_scheme(self):
        """Create the scheme to construct the global force vector.

        The scheme determines the initial force vector on all system dofs. 
        """
        self.main_model_part.ProcessInfo[StructuralMechanicsApplication.SYSTEM_DAMPING_RATIO] = self.harmonic_analysis_settings["system_damping"].GetDouble()
        self.main_model_part.ProcessInfo[StructuralMechanicsApplication.RAYLEIGH_ALPHA] = self.harmonic_analysis_settings["rayleigh_alpha"].GetDouble()
        self.main_model_part.ProcessInfo[StructuralMechanicsApplication.RAYLEIGH_BETA] = self.harmonic_analysis_settings["rayleigh_beta"].GetDouble()
        if self.settings["solution_type"].GetString() == "Dynamic":
            solution_scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
        else:
            raise Exception("Unsupported solution_type: " + self.settings["solution_type"])

        return solution_scheme

    def _create_linear_solver(self):
        """Create a dummy linear solver.
        
        This overrides the base class method and returns an empty linear solver as the harmonic
        analysis does not need a linear solver.
        """
        return KratosMultiphysics.LinearSolver()

    def _create_mechanical_solver(self):
        eigen_scheme = self.get_solution_scheme() # The scheme defines the matrices of the eigenvalue problem.
        builder_and_solver = self.get_builder_and_solver() # The eigensolver is created here.
        computing_model_part = self.GetComputingModelPart()

        return StructuralMechanicsApplication.HarmonicAnalysisStrategy(computing_model_part,
                                                                    eigen_scheme,
                                                                    builder_and_solver)
