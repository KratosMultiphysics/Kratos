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
        self.min_buffer_size = 2

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KM.MOMENTUM)
        self.main_model_part.AddNodalSolutionStepVariable(SW.ATMOSPHERIC_PRESSURE)

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.MOMENTUM_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.MOMENTUM_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.HEIGHT, self.main_model_part)

        KM.Logger.PrintInfo(self.__class__.__name__, "Shallow water solver DOFs added correctly.")

    def Initialize(self):
        super().Initialize()
        self.main_model_part.ProcessInfo.SetValue(SW.LUMPED_MASS_FACTOR, self.settings["lumped_mass_factor"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(SW.SHOCK_STABILIZATION_FACTOR, self.settings["shock_stabilization_factor"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(SW.GROUND_IRREGULARITY, self.settings["ground_irregularity"].GetDouble())

    def InitializeSolutionStep(self):
        SW.ShallowWaterUtilities().IdentifyWetDomain(self.main_model_part, KM.ACTIVE, self.main_model_part.ProcessInfo.GetValue(SW.DRY_HEIGHT))
        super().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        KM.VariableUtils().SetFlag(KM.ACTIVE, True, self.main_model_part.Nodes)
        KM.VariableUtils().SetFlag(KM.ACTIVE, True, self.main_model_part.Elements)
        dry_height = self.main_model_part.ProcessInfo.GetValue(SW.DRY_HEIGHT)
        SW.ShallowWaterUtilities().ComputeFreeSurfaceElevation(self.main_model_part)
        SW.ShallowWaterUtilities().ResetDryDomain(self.main_model_part, dry_height)
        SW.ComputeVelocityProcess(self.main_model_part, dry_height).Execute()
        self._CheckWaterLoss()

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KM.Parameters("""
        {
        "lumped_mass_factor"         : 1.0,
        "shock_stabilization_factor" : 0.001,
        "ground_irregularity"        : 0.0
        }
        """)
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

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
