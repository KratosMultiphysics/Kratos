# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

def CreateSolver(model, custom_settings):
    return HarmonicAnalysisSolver(model, custom_settings)

class HarmonicAnalysisSolver(MechanicalSolver):
    """The structural mechanics harmonic analysis solver.

    This class creates the mechanical solvers for the harmonic analysis.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[HarmonicAnalysisSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "scheme_type"   : "dynamic",
            "harmonic_analysis_settings" : {
                "use_effective_material_damping" : false
            }
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    #### Private functions ####

    def _CreateScheme(self):
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

    def _CreateLinearSolver(self):
        """Create a dummy linear solver.

        This overrides the base class method and returns an empty linear solver as the harmonic
        analysis does not need a linear solver.
        """
        return KratosMultiphysics.LinearSolver()

    def _CreateSolutionStrategy(self):
        eigen_scheme = self._GetScheme()
        builder_and_solver = self._GetBuilderAndSolver()
        computing_model_part = self.GetComputingModelPart()

        return StructuralMechanicsApplication.HarmonicAnalysisStrategy(computing_model_part,
                                                                    eigen_scheme,
                                                                    builder_and_solver,
                                                                    self.settings["harmonic_analysis_settings"]["use_effective_material_damping"].GetBool())
