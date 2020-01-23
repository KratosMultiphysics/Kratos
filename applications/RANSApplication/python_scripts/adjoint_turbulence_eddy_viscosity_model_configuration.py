import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS

from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
if CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    from KratosMultiphysics.FluidDynamicsApplication.adjoint_turbulence_model_solver import AdjointTurbulenceModelSolver
else:
    msg = "RANSApplication requires FluidDynamicsApplication which is not found."
    msg += " Please install/compile it and try again."
    raise Exception(msg)


class AdjointTurbulenceEddyViscosityModelConfiguration(
        AdjointTurbulenceModelSolver):
    def __init__(self, model, settings):
        self._validate_settings_in_baseclass = True  # To be removed eventually

        super(AdjointTurbulenceEddyViscosityModelConfiguration, self).__init__(
            model, settings)

        self.element_name = None
        self.condition_name = None

    def GetDefaultSettings(self):
        return Kratos.Parameters(r'''{
            "model_type"            : "",
            "model_settings"        : {},
            "echo_level"              : 0
        }''')

    def Initialize(self):
        super(AdjointTurbulenceEddyViscosityModelConfiguration,
              self).Initialize()

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                "Initialization successfull.")

    def AddVariables(self):
        # adding variables required by rans eddy viscosity models
        self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.DISTANCE)
        self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.KINEMATIC_VISCOSITY)
        self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.TURBULENT_VISCOSITY)
        self.fluid_model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_Y_PLUS)
        self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.RELAXED_ACCELERATION)

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                "Successfully added solution step variables.")

