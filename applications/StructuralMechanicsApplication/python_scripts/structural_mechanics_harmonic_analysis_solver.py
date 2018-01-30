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
                "use_effective_material_damping" : false
            }
        }
        """)
        self.validate_and_transfer_matching_settings(custom_settings, harmonic_analysis_settings)
        self.harmonic_analysis_settings = harmonic_analysis_settings["harmonic_analysis_settings"]
        # Validate the remaining settings in the base class.
        if not custom_settings.Has("scheme_type"): # Override defaults in the base class.
            custom_settings.AddEmptyValue("scheme_type")
            custom_settings["scheme_type"].SetString("dynamic")
        
        # Construct the base solver.
        super(HarmonicAnalysisSolver, self).__init__(main_model_part, custom_settings)
        print("::[HarmonicAnalysisSolver]:: Construction finished")

    #### Private functions ####

    def _create_solution_scheme(self):
        """Create the scheme to construct the global force vector.

        The scheme determines the initial force vector on all system dofs. 
        """
        if self.settings["scheme_type"].GetString() == "dynamic":
            solution_scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
        else:
            err_msg =  "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"dynamic\""
            raise Exception(err_msg)

        return solution_scheme

    def _create_linear_solver(self):
        """Create a dummy linear solver.
        
        This overrides the base class method and returns an empty linear solver as the harmonic
        analysis does not need a linear solver.
        """
        return KratosMultiphysics.LinearSolver()

    def _create_mechanical_solver(self):
        eigen_scheme = self.get_solution_scheme()
        builder_and_solver = self.get_builder_and_solver()
        computing_model_part = self.GetComputingModelPart()

        return StructuralMechanicsApplication.HarmonicAnalysisStrategy(computing_model_part,
                                                                    eigen_scheme,
                                                                    builder_and_solver,
                                                                    self.harmonic_analysis_settings["use_effective_material_damping"].GetBool())
