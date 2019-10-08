from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
#import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Import base class file
from KratosMultiphysics.RomApplication import ConvectionDiffusionROMsolver
import KratosMultiphysics.RomApplication as romapp

def CreateSolver(main_model_part, custom_settings):
    return HROMSolver(main_model_part, custom_settings)

class HROMSolver(ConvectionDiffusionROMsolver.ROMSolver):
    """The stationary class for Hyper Reduced Oorder Model (H-ROM) convection-diffusion solvers.

    Public member variables:
    stationary_settings -- settings for the implicit dynamic solvers.

    See convection_diffusion_base_solver.py for more information.
    """

    def __init__(self, main_model_part, custom_settings):
        # Set defaults and validate custom settings.
        self.stationary_settings = KratosMultiphysics.Parameters(r"""{}""")

        # Construct the base solver and validate the remaining settings in the base class
        super(HROMSolver, self).__init__(main_model_part, custom_settings)

        # Overwrite the base solver minimum buffer size
        if (self.settings["element_replace_settings"]["element_name"].GetString() == "EulerianConvDiff"):
            self.min_buffer_size = 2
        else:
            self.min_buffer_size = 1

        KratosMultiphysics.Logger.PrintInfo("::[HROMSolver]:: ", "Construction finished")


    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        rom_parameters=self.settings["rom_settings"]
        builder_and_solver = romapp.HROMBuilderAndSolver(linear_solver, rom_parameters)
        return builder_and_solver

    
    # def PrepareModelPart(self):
    #     if not self.is_restarted():
    #         # Check and prepare computing model part and import constitutive laws.
    #         self._execute_after_reading()

    #         throw_errors = False
    #         #KratosMultiphysics.TetrahedralMeshOrientationCheck(self.main_model_part, throw_errors).Execute()

    #         KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part,self._get_element_condition_replace_settings()).Execute()

    #         self._set_and_fill_buffer()

    #     if (self.settings["echo_level"].GetInt() > 0):
    #         KratosMultiphysics.Logger.PrintInfo(self.model)

    #     KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]::", "ModelPart prepared for Solver.")
