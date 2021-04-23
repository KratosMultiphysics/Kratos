# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_solver


def CreateSolver(model, custom_settings):
    return ConvectionDiffusionExplicitSolver(model, custom_settings)


class ConvectionDiffusionExplicitSolver(convection_diffusion_solver.ConvectionDiffusionSolver):
    """
    The explicit class for convection-diffusion solvers.
    See convection_diffusion_solver.py for more information.
    """

    def __init__(self, model, custom_settings):
        # Construct the base solver and validate the remaining settings in the base class
        super().__init__(model, custom_settings)

        # Overwrite the base solver minimum buffer size
        self.min_buffer_size = 2

        element_name = self.settings["element_replace_settings"]["element_name"].GetString()
        if self.settings["use_orthogonal_subscales"].GetBool() is True:
            oss_element_list = ["QSConvectionDiffusionExplicit", "DConvectionDiffusionExplicit"]
            if element_name in oss_element_list:
                self.main_model_part.ProcessInfo.SetValue(
                    KratosMultiphysics.OSS_SWITCH, 1
                )
            else:
                err_msg = (
                    "The selected element",
                    element_name,
                    "does not support OSS projection. Select QSConvectionDiffusionExplicit or DConvectionDiffusionExplicit instead.",
                )
                raise Exception(err_msg)
        else:
            if element_name in (
                "QSConvectionDiffusionExplicit",
                "DConvectionDiffusionExplicit",
            ):
                self.main_model_part.ProcessInfo.SetValue(
                    KratosMultiphysics.OSS_SWITCH, 0
                )

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KratosMultiphysics.Parameters(
        """
        {
            "time_integration_method" : "explicit",
            "use_orthogonal_subscales" : false,
            "dynamic_tau": 1.0
        }
        """
        )

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def Initialize(self):
        super().Initialize()
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DYNAMIC_TAU] = self.settings["dynamic_tau"].GetDouble()

    #### Private functions ####

    def _create_builder_and_solver(self):
        builder_and_solver = KratosMultiphysics.ExplicitBuilder()
        return builder_and_solver

    def _create_convection_diffusion_solution_strategy(self):
        convection_diffusion_solution_strategy = self._create_runge_kutta_4_strategy()
        return convection_diffusion_solution_strategy

    def _create_runge_kutta_4_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        explicit_builder_and_solver = self.get_builder_and_solver()
        rebuild_level = 0
        return ConvectionDiffusionApplication.ExplicitSolvingStrategyRungeKutta4ConvectionDiffusion(
            computing_model_part,
            explicit_builder_and_solver,
            self.settings["move_mesh_flag"].GetBool(),
            rebuild_level,
        )
