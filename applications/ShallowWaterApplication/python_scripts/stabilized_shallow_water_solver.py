# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return StabilizedShallowWaterSolver(model, custom_settings)

class StabilizedShallowWaterSolver(ShallowWaterBaseSolver):
    def __init__(self, model, settings):
        super().__init__(model, settings)

        # Set the element and condition names for the replace settings
        self.element_name = "ShallowWater"
        self.condition_name = "LineCondition"
        self.min_buffer_size = self.settings["time_integration_order"].GetInt() + 1

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KM.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(SW.VERTICAL_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(SW.ATMOSPHERIC_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KM.MESH_ACCELERATION)

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.MOMENTUM_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.MOMENTUM_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.HEIGHT, self.main_model_part)
        KM.Logger.PrintInfo(self.__class__.__name__, "Shallow water solver DOFs added correctly.")

    def Initialize(self):
        super().Initialize()
        self.main_model_part.ProcessInfo.SetValue(SW.RELATIVE_DRY_HEIGHT, self.settings["relative_dry_height"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KM.STABILIZATION_FACTOR, self.settings["stabilization_factor"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(SW.SHOCK_STABILIZATION_FACTOR, self.settings["shock_stabilization_factor"].GetDouble())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        SW.ShallowWaterUtilities().ComputeFreeSurfaceElevation(self.main_model_part)
        SW.ShallowWaterUtilities().ComputeVelocity(self.main_model_part, True)
        self._CheckWaterLoss()

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KM.Parameters("""
        {
        "time_integration_order"     : 2,
        "relative_dry_height"        : 0.1,
        "stabilization_factor"       : 0.005,
        "shock_stabilization_factor" : 0.001,
        "add_flux_correction"        : false
        }
        """)
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def _CreateScheme(self):
        if self.settings["add_flux_correction"].GetBool():
            time_scheme = SW.FluxCorrectedShallowWaterScheme(self.settings["time_integration_order"].GetInt())
            if self.settings["shock_stabilization_factor"].GetDouble() > 0.0:
                KM.Logger.PrintWarning(self.__class__.__name__, "Detected shock stabilization with flux correction, please, disable on of them.")
        else:
            time_scheme = SW.ShallowWaterResidualBasedBDFScheme(self.settings["time_integration_order"].GetInt())
        return time_scheme

    def _InitializeWaterLoss(self):
        self.initial_water = KM.VariableUtils().SumHistoricalNodeScalarVariable(SW.HEIGHT, self.main_model_part,0)
        self.initial_water /= self.main_model_part.NumberOfNodes()

    def _CheckWaterLoss(self):
        if not hasattr(self, 'initial_water'):
            self._InitializeWaterLoss()
        total_water = KM.VariableUtils().SumHistoricalNodeScalarVariable(SW.HEIGHT, self.main_model_part,0)
        total_water /= self.main_model_part.NumberOfNodes()
        water_loss = (total_water - self.initial_water) / self.initial_water
        if abs(water_loss) > 1e-3 and self.echo_level > 1:
            msg = "Water loss : {} %"
            KM.Logger.PrintWarning(self.__class__.__name__, msg.format(water_loss*100))
