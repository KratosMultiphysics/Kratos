from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Import base class file
import convection_diffusion_stationary_solver

def CreateSolver(main_model_part, custom_settings):
    return ConvectionDiffusionStationarySolverTemporary(main_model_part, custom_settings)

class ConvectionDiffusionStationarySolverTemporary(convection_diffusion_stationary_solver.ConvectionDiffusionStationarySolver):
    """Derived class for convection-diffusion solvers.

    Public member variables:
    stationary_settings -- settings for the implicit dynamic solvers.

    See convection_diffusion_base_solver.py for more information.
    """

    def PrepareModelPart(self):
        if not self.is_restarted():
            # Check and prepare computing model part and import constitutive laws.
            self._execute_after_reading()

            # throw_errors = False
            # KratosMultiphysics.TetrahedralMeshOrientationCheck(self.main_model_part, throw_errors).Execute()

            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part,self._get_element_condition_replace_settings()).Execute()

            self._set_and_fill_buffer()

        if (self.settings["echo_level"].GetInt() > 0):
            self.print_on_rank_zero(self.model)

        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]::", "ModelPart prepared for Solver.")

