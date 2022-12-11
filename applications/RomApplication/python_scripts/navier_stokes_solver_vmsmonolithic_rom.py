# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_solver_vmsmonolithic import NavierStokesSolverMonolithic
import KratosMultiphysics.RomApplication as romapp

def CreateSolver(model, custom_settings):
    return ROMSolver(model, custom_settings)

class ROMSolver(NavierStokesSolverMonolithic):

    def __init__(self, model, custom_settings):
        KratosMultiphysics.Logger.PrintWarning('\x1b[1;31m[DEPRECATED CLASS] \x1b[0m',"\'navier_stokes_solver_vmsmonolithic_rom\'", "class is deprecated. Use the generic\'RomSolver\' one instead.")

        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[ROMSolver]:: ", "Construction finished")

    #### Private functions ####
    @classmethod
    def GetDefaultParameters(cls):
        import json
        with open("ProblemFiles/RomParameters.json", 'r') as file:
            f = json.load(file)
        with open('ProblemFiles/RomParametersHeader.json', 'w') as file:
            json.dump({"rom_settings":f["rom_settings"]},file, indent=2)
        with open("ProblemFiles/RomParametersHeader.json", 'r') as file:
            default_settings = KratosMultiphysics.Parameters(file.read())
        default_settings.AddMissingParameters(super(ROMSolver,cls).GetDefaultParameters())
        return default_settings

    def _CreateBuilderAndSolver(self):
        linear_solver = self._GetLinearSolver()
        rom_parameters=self.settings["rom_settings"]
        builder_and_solver = romapp.ROMBuilderAndSolver(linear_solver, rom_parameters)
        return builder_and_solver
