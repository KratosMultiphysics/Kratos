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
                "assembling_strategy" : "global",
                "rom_settings": {
                    "nodal_unknowns": [],
                    "number_of_rom_dofs": 0,
                    "rom_bns_settings": {}
                }
            }""")
            default_settings.AddMissingParameters(super().GetDefaultParameters())
            return default_settings

        def _CreateBuilderAndSolver(self):
            linear_solver = self._GetLinearSolver()
            rom_parameters, solving_strategy = self._ValidateAndReturnRomParameters()
            available_solving_strategies = {
                "elemental_galerkin": KratosROM.ROMBuilderAndSolver,
                "global_galerkin": KratosROM.GlobalROMBuilderAndSolver,
                "lspg": KratosROM.LeastSquaresPetrovGalerkinROMBuilderAndSolver,
                "elemental_petrov_galerkin": KratosROM.PetrovGalerkinROMBuilderAndSolver,
                "global_petrov_galerkin": KratosROM.GlobalPetrovGalerkinROMBuilderAndSolver
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
            assembling_strategy = self.settings["assembling_strategy"].GetString()

            # Check for the 'elemental' assembling strategy in case of 'lspg' projection strategy
            if projection_strategy == "lspg" and assembling_strategy == "elemental":
                warn_msg = "Elemental LSPG is not yet available. Using default global assembling strategy instead."
                KratosMultiphysics.Logger.PrintWarning("::[ROMSolver]:: ", warn_msg)

            # For now, only Galerkin and Petrov-Galerkin projections have the elemental or global approach option
            if projection_strategy in ("galerkin", "petrov_galerkin"): #TODO: Possibility of doing elemental lspg
                available_assembling_strategies = {
                    "global",
                    "elemental"
                }
                if assembling_strategy in available_assembling_strategies:
                    # Add the assembling strategy prefix to the projection strategy
                    projection_strategy = f"{assembling_strategy}_{projection_strategy}"
                else:
                    err_msg = f"'Assembling_strategy': '{assembling_strategy}' is not available. Please select one of the following: {list(available_assembling_strategies)}."
                    raise ValueError(err_msg)

            self._AssignMissingInnerRomParameters(projection_strategy)

            # Return the validated ROM parameters
            return self.settings["rom_settings"], projection_strategy

        def _AssignMissingInnerRomParameters(self, projection_strategy):
            monotonicity_preserving = self.settings["rom_settings"]["rom_bns_settings"]["monotonicity_preserving"].GetBool() if self.settings["rom_settings"]["rom_bns_settings"].Has("monotonicity_preserving") else False
            if projection_strategy in ("global_galerkin", "lspg", "global_petrov_galerkin"):
                self.settings["rom_settings"]["rom_bns_settings"].AddBool("monotonicity_preserving", monotonicity_preserving)

    return ROMSolver(model, custom_settings)