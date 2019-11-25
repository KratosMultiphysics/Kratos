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

    See structural_mechanics_solver.py for more information.
    """

    def __init__(self, main_model_part, custom_settings):
        # Set defaults and validate custom settings.
        #self.stationary_settings = KratosMultiphysics.Parameters(r"""{}""")

        # Construct the base solver and validate the remaining settings in the base class
        super(ROMSolver, self).__init__(main_model_part, custom_settings)

        # # Overwrite the base solver minimum buffer size
        # if (self.settings["element_replace_settings"]["element_name"].GetString() == "EulerianConvDiff"):
        #     self.min_buffer_size = 2
        # else:
        #     self.min_buffer_size = 1

        KratosMultiphysics.Logger.PrintInfo("::[ROMSolver]:: ", "Construction finished")

    #### Private functions ####
    @classmethod
    def GetDefaultSettings(cls):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "rom_settings": {
            "nodal_unknowns": [ "DISPLACEMENT", "DISPLACEMENT_Y"],
            "number_of_rom_dofs": 3
            }
        }
        """)
        default_settings.AddMissingParameters(super(ROMSolver,cls).GetDefaultSettings())
        return default_settings

    def AddVariables(self):
        super(ROMSolver, self).AddVariables() #Adding nodal area variable
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)


    # def _create_solution_scheme(self):
    #     #Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
    #     #self.GetComputingModelPart().ProcessInfo[ConvectionDiffusionApplication.THETA] = 1.0
    #     #self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DYNAMIC_TAU] = 0.0
    #     structural_mechanics_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
    #     return structural_mechanics_scheme


    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        rom_parameters=self.settings["rom_settings"]
        builder_and_solver = romapp.ROMBuilderAndSolver(linear_solver, rom_parameters)
        return builder_and_solver
