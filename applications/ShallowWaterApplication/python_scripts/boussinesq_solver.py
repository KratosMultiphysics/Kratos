# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return BoussinesqSolver(model, custom_settings)

class BoussinesqSolver(ShallowWaterBaseSolver):
    def __init__(self, model, settings):
        super().__init__(model, settings)
        self.element_name = "BoussinesqElement"
        self.condition_name = "BoussinesqCondition"
        self.min_buffer_size = 4

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.VELOCITY_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.VELOCITY_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.FREE_SURFACE_ELEVATION, self.main_model_part)
        KM.Logger.PrintInfo(self.__class__.__name__, "Boussinesq equations DOFs added correctly.")

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KM.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(SW.VERTICAL_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KM.VELOCITY_LAPLACIAN)
        self.main_model_part.AddNodalSolutionStepVariable(KM.RESIDUAL_VECTOR) # This is used by the predictor. TODO: replace with reaction
        self.main_model_part.AddNodalSolutionStepVariable(KM.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(SW.FIRST_DERIVATIVE_WEIGHTS)
        self.main_model_part.AddNodalSolutionStepVariable(SW.SECOND_DERIVATIVE_WEIGHTS)

    def AdvanceInTime(self, current_time):
        current_time = super().AdvanceInTime(current_time)
        if self._TimeBufferIsInitialized():
            current_time_step = self.GetComputingModelPart().ProcessInfo.GetValue(KM.DELTA_TIME)
            previous_time_step = self.GetComputingModelPart().ProcessInfo.GetPreviousTimeStepInfo().GetValue(KM.DELTA_TIME)
            if current_time_step - previous_time_step > 1e-10:
                KM.Logger.PrintWarning(self.__class__.__name__, "The Adams Moulton scheme requires a constant time step.")
        return current_time

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        SW.ShallowWaterUtilities().ComputeHeightFromFreeSurface(self.main_model_part)

    def _SetProcessInfo(self):
        super()._SetProcessInfo()
        self.main_model_part.ProcessInfo.SetValue(SW.AMPLITUDE, self.settings["amplitude"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(SW.WAVELENGTH, self.settings["wavelength"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KM.STABILIZATION_FACTOR, self.settings["stabilization_factor"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(SW.RELATIVE_DRY_HEIGHT, self.settings["relative_dry_height"].GetDouble())

    def _CreateScheme(self):
        return SW.ResidualBasedAdamsMoultonScheme()

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KM.Parameters("""{
            "relative_dry_height"        : 0.1,
            "stabilization_factor"       : 0.01,
            "wavelength"                 : 10,
            "amplitude"                  : 0.0
        }""")
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings
