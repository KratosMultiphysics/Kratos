# Importing the Kratos Library
from typing import List
import KratosMultiphysics

# Import applications
import KratosMultiphysics.RomApplication as KratosROM

def CreateSolver(cls, model, custom_settings):
    class ROMSolver(cls):
        """ROM solver generic main class.
        This class serves as a generic class to make a standard Kratos solver ROM-compatible.
        It extends the default parameters to include the \'rom_settings\' and overrides the
        creation of the builder and solver to use the ROM one.
        """

        def __init__(self, model, custom_settings):
            super().__init__(model, custom_settings)
            KratosMultiphysics.Logger.PrintInfo("::[ROMSolver]:: ", "Construction finished")

        @classmethod
        def GetDefaultParameters(cls):
            default_settings = KratosMultiphysics.Parameters("""{
                "projection_strategy" : "galerkin",
                "rom_settings": {
                    "nodal_unknowns": [],
                    "number_of_rom_dofs": 0
                }
            }""")
            default_settings.AddMissingParameters(super().GetDefaultParameters())
            return default_settings

        def _CreateBuilderAndSolver(self):
            linear_solver = self._GetLinearSolver()
            rom_parameters, solving_strategy = self._ValidateAndReturnRomParameters()
            available_solving_strategies = {
                "galerkin": KratosROM.ROMBuilderAndSolver, 
                "lspg": KratosROM.LeastSquaresPetrovGalerkinROMBuilderAndSolver,
                "petrov_galerkin": KratosROM.PetrovGalerkinROMBuilderAndSolver
            }
            if solving_strategy in available_solving_strategies:
                return available_solving_strategies[solving_strategy](linear_solver, rom_parameters)
            else:
                err_msg = f"'Solving_strategy': '{solving_strategy}' is not available. Please select one of the following: {list(available_solving_strategies.keys())}."
                raise ValueError(err_msg)

        def _ValidateAndReturnRomParameters(self):
            # Check that the number of ROM DOFs has been provided
            n_rom_dofs = self.settings["rom_settings"]["number_of_rom_dofs"].GetInt()
            if not n_rom_dofs > 0:
                err_msg = "\'number_of_rom_dofs\' in \'rom_settings\' is {}. Please set a larger than zero value.".format(n_rom_dofs)
                raise Exception(err_msg)

            # Check if the nodal unknowns have been provided by the user
            # If not, take the DOFs list from the base solver
            nodal_unknowns = self.settings["rom_settings"]["nodal_unknowns"].GetStringArray()
            if len(nodal_unknowns) == 0:
                solver_dofs_list = self.GetDofsList()
                if not len(solver_dofs_list) == 0:
                    self.settings["rom_settings"]["nodal_unknowns"].SetStringArray(solver_dofs_list)
                else:
                    err_msg = "\'nodal_unknowns\' in \'rom_settings\' is not provided and there is a not-valid implementation in base solver."
                    err_msg += " Please manually set \'nodal_unknowns\' in \'rom_settings\'."
                    raise Exception(err_msg)
            projection_strategy = self.settings["projection_strategy"].GetString()

            # Return the validated ROM parameters
            return self.settings["rom_settings"], projection_strategy

    return ROMSolver(model, custom_settings)