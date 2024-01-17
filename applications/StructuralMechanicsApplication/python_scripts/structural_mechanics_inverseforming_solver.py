# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

def CreateSolver(model, custom_settings):
    return InverseFormingMechanicalSolver(model, custom_settings)

class InverseFormingMechanicalSolver(MechanicalSolver):
    """The structural mechanics inverse forming solver.

    This class creates the mechanical solver for inverseforming.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[InverseFormingMechanicalSolver]:: ", "Construction finished")
    
    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "printing_format"             : "none"
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults
    
    def _CreateScheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
    
    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._GetScheme()
        mechanical_convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        return StructuralMechanicsApplication.InverseFormingStrategy(
                                                                    computing_model_part,
                                                                    mechanical_scheme,
                                                                    mechanical_convergence_criterion,
                                                                    builder_and_solver,
                                                                    self.settings["max_iteration"].GetInt(),
                                                                    self.settings["compute_reactions"].GetBool(),
                                                                    self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                    self.settings["move_mesh_flag"].GetBool())