from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_base_solver
from StructuralMechanicsROMsolver import ROMSolver
import KratosMultiphysics.RomApplication as romapp

def CreateSolver(main_model_part, custom_settings):
    return ROMSolver(main_model_part, custom_settings)

class HROMSolver(ROMSolver):
    """The stationary class for H-ROM structural mechanics solvers.

    Public member variables:
    stationary_settings -- settings for the implicit dynamic solvers.

    See structural_mechanics_solver.py for more information.
    """

    def __init__(self, main_model_part, custom_settings):
        super(HROMSolver, self).__init__(main_model_part, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[HROMSolver]:: ", "Construction finished")


    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        rom_parameters=self.settings["rom_settings"]
        builder_and_solver = romapp.HROMBuilderAndSolver(linear_solver, rom_parameters)
        return builder_and_solver
