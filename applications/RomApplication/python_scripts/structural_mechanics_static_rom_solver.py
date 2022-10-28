
# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_static_solver import StaticMechanicalSolver
import KratosMultiphysics.RomApplication as romapp

def CreateSolver(main_model_part, custom_settings):
    return ROMSolver(main_model_part, custom_settings)

class ROMSolver(StaticMechanicalSolver):
    """The stationary class for ROM structural mechanics solvers.

    See structural_mechanics_static_solver.py for more information.
    """

    def __init__(self, main_model_part, custom_settings):
        KratosMultiphysics.Logger.PrintWarning('\x1b[1;31m[DEPRECATED CLASS] \x1b[0m',"\'stuctural_mechanics_static_rom_solver\'", "class is deprecated. Use the generic\'RomSolver\' one instead.")

        super(ROMSolver, self).__init__(main_model_part, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[ROMSolver]:: ", "Construction finished")

    #### Private functions ####
    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "rom_settings": {
            "nodal_unknowns": [ "DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z"],
            "number_of_rom_dofs": 3
            }
        }
        """)
        default_settings.AddMissingParameters(super(ROMSolver,cls).GetDefaultParameters())
        return default_settings

    def AddVariables(self):
        super(ROMSolver, self).AddVariables() #Adding nodal area variable
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

    def _CreateBuilderAndSolver(self):
        linear_solver = self._GetLinearSolver()
        rom_parameters=self.settings["rom_settings"]
        builder_and_solver = romapp.ROMBuilderAndSolver(linear_solver, rom_parameters)
        return builder_and_solver
