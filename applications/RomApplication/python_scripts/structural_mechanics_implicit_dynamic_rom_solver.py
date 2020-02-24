from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_implicit_dynamic_solver import ImplicitMechanicalSolver
import KratosMultiphysics.RomApplication as romapp

def CreateSolver(main_model_part, custom_settings):
    return ROMSolver(main_model_part, custom_settings)

class ROMSolver(ImplicitMechanicalSolver):
    """The stationary class for ROM structural mechanics solvers.

    Public member variables:
    stationary_settings -- settings for the implicit dynamic solvers.

    See structural_mechanics_implicit_dynamic_solver.py for more information.
    """

    def __init__(self, main_model_part, custom_settings):
        super(ROMSolver, self).__init__(main_model_part, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[ROMSolver]:: ", "Construction finished")

    #### Private functions ####
    @classmethod
    def GetDefaultSettings(cls):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "rom_settings": {
            "nodal_unknowns": [ "DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z"],
            "number_of_rom_dofs": 10
            }
        }
        """)
        default_settings.AddMissingParameters(super(ROMSolver,cls).GetDefaultSettings())
        return default_settings

    def AddVariables(self):
        super(ROMSolver, self).AddVariables() #Adding nodal area variable
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        rom_parameters=self.settings["rom_settings"]
        builder_and_solver = romapp.ROMBuilderAndSolver(linear_solver, rom_parameters)
        return builder_and_solver
