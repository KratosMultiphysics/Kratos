# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

def CreateSolver(model, custom_settings):
    return StaticMechanicalSolver(model, custom_settings)

class StaticMechanicalSolver(MechanicalSolver):
    """The structural mechanics static solver.

    This class creates the mechanical solvers for static analysis.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[StaticMechanicalSolver]:: ", "Construction finished")

    def _CreateScheme(self):
        scheme_settings = KratosMultiphysics.Parameters("""{
            "projection_variables_list" : []
        }""")
        if self.settings["use_orthogonal_subscales"].GetBool():
            if self.settings["volumetric_strain_dofs"].GetBool():
                scheme_settings["projection_variables_list"].SetStringArray(["VOLUMETRIC_STRAIN_PROJECTION","DISPLACEMENT_PROJECTION"])
        return KratosStructural.StructuralMechanicsStaticScheme(scheme_settings)
