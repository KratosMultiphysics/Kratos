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

        self._SetUpFormulation()
        self.min_buffer_size = self.settings["time_integration_order"].GetInt() + 1

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KM.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(SW.VERTICAL_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(SW.ATMOSPHERIC_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KM.MESH_ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(SW.WIND)
        self.main_model_part.AddNodalSolutionStepVariable(KM.DISTANCE)

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
        self.main_model_part.ProcessInfo.SetValue(KM.DENSITY_AIR, 1e0)
        self.main_model_part.ProcessInfo.SetValue(KM.DENSITY, 1e3)
        self.main_model_part.ProcessInfo.SetValue(SW.INTEGRATE_BY_PARTS, False)
        if self.compute_neighbors:
            dimension = self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE]
            avg_neigh = 2 * dimension
            KM.FindElementalNeighboursProcess(self.main_model_part, dimension, avg_neigh).Execute()

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
        "stabilization_factor"       : 0.01,
        "shock_stabilization_factor" : 1.0,
        "shock_capturing_type"       : "residual_viscosity"
        }
        """)
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def _CreateScheme(self):
        if self.add_flux_correction:
            time_scheme = SW.FluxCorrectedShallowWaterScheme(self.settings["time_integration_order"].GetInt())
            if self.settings["shock_stabilization_factor"].GetDouble() > 0.0:
                self.settings["shock_stabilization_factor"].SetDouble(0.0)
                KM.Logger.PrintWarning(self.__class__.__name__, "Detected shock stabilization with flux correction. The shock stabilization factor will be set to 0.")
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

    def _SetUpFormulation(self):
        shock_capturing_type = self.settings["shock_capturing_type"].GetString()
        if  shock_capturing_type == "residual_viscosity":
            self.element_name = "ShallowWater"
            self.condition_name = "LineCondition"
            self.compute_neighbors = False
            self.add_flux_correction = False
        elif shock_capturing_type == "flux_correction":
            self.element_name = "ShallowWater"
            self.condition_name = "LineCondition"
            self.compute_neighbors = False
            self.add_flux_correction = True
        elif shock_capturing_type == "gradient_jump":
            self.element_name = "ConservativeElement"
            self.condition_name = "ConservativeCondition"
            self.compute_neighbors = True
            self.add_flux_correction = False
        else:
            msg  = "StabilizedShallowWaterSolver._SetUpFormulation:\n"
            msg += "The specified 'shock_capturing_type' : '{}' is not available.\n".format(shock_capturing_type)
            msg += "The possible options are:\n"
            msg += "\t- 'residual_viscosity'\n"
            msg += "\t- 'flux_correction'\n"
            msg += "\t- 'gradient_jump'\n"
            raise Exception(msg)
