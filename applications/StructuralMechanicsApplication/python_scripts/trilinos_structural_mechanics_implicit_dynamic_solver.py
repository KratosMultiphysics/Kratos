# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.TrilinosApplication as TrilinosApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.trilinos_structural_mechanics_solver import TrilinosMechanicalSolver

def CreateSolver(model, custom_settings):
    return TrilinosImplicitMechanicalSolver(model, custom_settings)

class TrilinosImplicitMechanicalSolver(TrilinosMechanicalSolver):
    """The trilinos structural mechanics implicit dynamic solver.

    For more information see:
    structural_mechanics_solver.py
    trilinos_structural_mechanics_solver.py
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super(TrilinosImplicitMechanicalSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosImplicitMechanicalSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "time_integration_method" : "implicit",
            "scheme_type"             : "bossak",
            "damp_factor_m"           :-0.3,
            "rayleigh_alpha"          : 0.0,
            "rayleigh_beta"           : 0.0
        }""")
        this_defaults.AddMissingParameters(super(TrilinosImplicitMechanicalSolver, cls).GetDefaultSettings())
        return this_defaults

    def AddVariables(self):
        super(TrilinosImplicitMechanicalSolver, self).AddVariables()
        self._add_dynamic_variables()
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosImplicitMechanicalSolver]:: Variables ADDED")

    def AddDofs(self):
        super(TrilinosImplicitMechanicalSolver, self).AddDofs()
        self._add_dynamic_dofs()
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosImplicitMechanicalSolver]:: DOF's ADDED")

    #### Private functions ####

    def _create_solution_scheme(self):
        scheme_type = self.settings["scheme_type"].GetString()
        self.main_model_part.ProcessInfo[StructuralMechanicsApplication.RAYLEIGH_ALPHA] = self.settings["rayleigh_alpha"].GetDouble()
        self.main_model_part.ProcessInfo[StructuralMechanicsApplication.RAYLEIGH_BETA] = self.settings["rayleigh_beta"].GetDouble()
        if (scheme_type == "newmark"):
            damp_factor_m = 0.0
        elif (scheme_type == "bossak"):
            damp_factor_m = self.settings["damp_factor_m"].GetDouble()
        else:
            err_msg =  "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"newmark\", \"bossak\""
            raise Exception(err_msg)
        mechanical_scheme = TrilinosApplication.TrilinosResidualBasedBossakDisplacementScheme(damp_factor_m)
        return mechanical_scheme
