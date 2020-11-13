# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return FreeSurfaceShallowWaterSolver(model, custom_settings)

class FreeSurfaceShallowWaterSolver(ShallowWaterBaseSolver):
    def __init__(self, model, settings):
        super().__init__(model, settings)

        # Set the element and condition names for the replace settings
        self.element_name = "SWE"
        self.condition_name = "LineCondition"
        self.min_buffer_size = 2
        self.advection_epsilon = self.settings["advection_epsilon"].GetDouble()

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(SW.EQUIVALENT_MANNING)
        self.main_model_part.AddNodalSolutionStepVariable(KM.POROSITY)

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.MOMENTUM_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.MOMENTUM_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.FREE_SURFACE_ELEVATION, self.main_model_part)

        KM.Logger.PrintInfo(self.__class__.__name__, "Shallow water solver DOFs added correctly.")

    def InitializeSolutionStep(self):
        self._GetSolutionStrategy().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        epsilon = max(self.advection_epsilon, self.GetComputingModelPart().ProcessInfo[SW.DRY_HEIGHT])
        SW.ShallowWaterUtilities().ComputeHeightFromFreeSurface(self.GetComputingModelPart())
        SW.ComputeVelocityProcess(self.GetComputingModelPart(), epsilon).Execute()
        SW.ShallowWaterUtilities().ComputeAccelerations(self.GetComputingModelPart())

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KM.Parameters("""
        {
            "time_scale"            : "seconds",
            "water_height_scale"    : "meters",
            "advection_epsilon"     : 1.0e-2,
            "permeability"          : 1.0e-4,
            "dry_discharge_penalty" : 1.0e+2,
            "wetting_drying_model"  : {
                "model_name"    : "negative_height",
                "beta"          : 1e4
            }
        }""")
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def PrepareModelPart(self):
        super().PrepareModelPart()
        permeability = self.settings["permeability"].GetDouble()
        discharge_penalty = self.settings["dry_discharge_penalty"].GetDouble()
        time_scale = self.settings["time_scale"].GetString()
        water_height_scale = self.settings["water_height_scale"].GetString()

        # Checks
        if permeability == 0.0:
            KM.Logger.PrintWarning(self.__class__.__name__, "Detected permeability == 0.0")
        if discharge_penalty == 0.0:
            KM.Logger.PrintWarning(self.__class__.__name__, "Detected dry_discharge_penalty == 0.0")

        # Time unit converter
        if   time_scale == "seconds":
            time_unit_converter =     1
        elif time_scale == "minutes":
            time_unit_converter =    60
        elif time_scale == "hours":
            time_unit_converter =  3600
        elif time_scale == "days":
            time_unit_converter = 86400
        else:
            raise Exception("unknown time scale")

        # Water height unit converter
        if   water_height_scale == "meters":
            water_height_unit_converter = 1.0
        elif water_height_scale == "millimeters":
            water_height_unit_converter = 0.001
        else:
            raise Exception("unknown water height scale")

        self.main_model_part.ProcessInfo.SetValue(SW.PERMEABILITY, permeability)
        self.main_model_part.ProcessInfo.SetValue(SW.DRY_DISCHARGE_PENALTY, discharge_penalty)
        self.main_model_part.ProcessInfo.SetValue(SW.TIME_UNIT_CONVERTER, time_unit_converter)
        self.main_model_part.ProcessInfo.SetValue(SW.WATER_HEIGHT_UNIT_CONVERTER, water_height_unit_converter)

    def _GetWettingModel(self):
        if not hasattr(self, "_dry_wet_model"):
            self._wetting_model = self._CreateWettingModel()
        return self._wetting_model

    def _CreateWettingModel(self):
        if self.settings["wetting_drying_model"].Has("model_name"):
            if self.settings["wetting_drying_model"]["model_name"].GetString() == "rough_porous_layer":
                return SW.RoughPorousLayerWettingModel(self.GetComputingModelPart(), self.settings["wetting_drying_model"])
            if self.settings["wetting_drying_model"]["model_name"].GetString() == "negative_height":
                return SW.NegativeHeightWettingModel(self.GetComputingModelPart(), self.settings["wetting_drying_model"])
            else:
                msg = "Requested wetting drying model: " + self.settings["wetting_drying_model"]["model_name"].GetString() +"\n"
                msg += "Available options are:\n"
                msg += "\t\"rough_porous_layer\"\n"
                msg += "\t\"negative_height\"\n"
                raise Exception(msg)
        else:
            return None

    def _CreateScheme(self):
        wetting_model = self._GetWettingModel()
        time_scheme = SW.ResidualBasedIncrementalUpdateWettingScheme(wetting_model)
        return time_scheme
