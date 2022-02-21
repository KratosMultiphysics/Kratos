# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.RomApplication as RomApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

# Other imports

def CreateSolver(model, custom_settings):
    return ModalDerivativeSolver(model, custom_settings)

class ModalDerivativeSolver(MechanicalSolver):
    """The structural mechanics modal derivative solver.

    This class creates the modal derivative solver for computing modal derivatives.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        """Initializes the solver."""
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[ModalDerivativeSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "derivative_type"               : "static",
            "finite_difference_type"        : "forward",
            "finite_difference_step_size"   : 1e-3,
            "mass_orthonormalize"           : true,
            "rom_settings"                  :
            {
                "rom_parameters_filename"       : "RomParameters.json",
                "nodal_unknowns"                : []
            }
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def _CreateScheme(self):
        return RomApplication.ModalDerivativeScheme(self.settings)

    def _CreateBuilderAndSolver(self):

        if not self.settings["builder_and_solver_settings"]["use_block_builder"].GetBool():
            err_msg  = 'Modal derivative analysis is only available with block builder and solver!'
            raise Exception(err_msg)

        return super()._CreateBuilderAndSolver()

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        modal_derivative_scheme = self._GetScheme()
        builder_and_solver = self._GetBuilderAndSolver()

        return RomApplication.ModalDerivativeStrategy(computing_model_part,
                                                          modal_derivative_scheme,
                                                          builder_and_solver,
                                                          self.settings)
    
