# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.wave_solver import WaveSolver

def CreateSolver(model, custom_settings):
    return EmbeddedPrimitiveSolver(model, custom_settings)


class EmbeddedPrimitiveSolver(WaveSolver):
    def __init__(self, model, settings):
        super().__init__(model, settings)
        self.element_name, self.condition_name, self.min_buffer_size = self._GetFormulationSettings()

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KM.DISTANCE)

    def Initialize(self):
        super().Initialize()

        # Set the distance modification process
        self.GetDistanceModificationProcess().ExecuteInitialize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        # Perform levelset convection
        #TODO: To be implemented

        # Do the distance correction and set the fixity in the "internal" nodes
        self.GetDistanceModificationProcess().ExecuteInitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        # Restore the fixity to its original status
        self.GetDistanceModificationProcess().ExecuteFinalizeSolutionStep()

    def _GetFormulationSettings(self):
        scheme = self.settings["time_integration_scheme"].GetString()
        order = self.settings["time_integration_order"].GetInt()
        if scheme == "bdf":
            element_name = "EmbeddedPrimitiveElement"
            condition_name = "WaveCondition"
            buffer_size = order + 1
        else:
            raise Exception('The possible "time_integration_scheme" are "bdf" and "crank_nicolson"')
        return element_name, condition_name, buffer_size

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KM.Parameters("""{
            "distance_modification_settings" : {}
        }""")
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    @classmethod
    def __GetDistanceModificationDefaultSettings(self):
        return KM.Parameters(r'''{
            "model_part_name": "",
            "distance_threshold": 1e-3,
            "continuous_distance": true,
            "check_at_each_time_step": true,
            "avoid_almost_empty_elements": true,
            "deactivate_full_negative_elements": true,
            "full_negative_elements_fixed_variables_list" : ["HEIGHT","VELOCITY"]
        }''')

    def __CreateDistanceModificationProcess(self):
        # Set the distance modification settings according to the level set type
        # Note that the distance modification process is applied to the volume model part
        distance_modification_settings = self.settings["distance_modification_settings"]
        distance_modification_settings.ValidateAndAssignDefaults(
            self.__GetDistanceModificationDefaultSettings())
        volume_part_name = self.settings["model_part_name"].GetString()
        distance_modification_settings["model_part_name"].SetString(volume_part_name)
        return KratosFluid.DistanceModificationProcess(self.model, distance_modification_settings)

    def GetDistanceModificationProcess(self):
        if not hasattr(self, '_distance_modification_process'):
            self._distance_modification_process = self.__CreateDistanceModificationProcess()
        return self._distance_modification_process
