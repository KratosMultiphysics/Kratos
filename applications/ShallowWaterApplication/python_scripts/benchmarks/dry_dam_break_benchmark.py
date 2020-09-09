# Importing the Kratos Library
import KratosMultiphysics as KM

from KratosMultiphysics.ShallowWaterApplication.benchmarks.base_benchmark_process import BaseBenchmarkProcess

# Other imports
import numpy as np

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DryDamBreakBenchmark(model, settings["Parameters"])

class DryDamBreakBenchmark(BaseBenchmarkProcess):

    def __init__(self, model, settings):
        # The base class sets the model_part, variables and benchmark_settings
        super(DryDamBreakBenchmark, self).__init__(model, settings)

        benchmark_default_settings = KM.Parameters("""
            {
                "dam_position"  : 5.0,
                "left_height"   : 2.0
            }
            """
            )
        self.benchmark_settings.ValidateAndAssignDefaults(benchmark_default_settings)

        self.dam = self.benchmark_settings["dam_position"].GetDouble()
        self.hl = self.benchmark_settings["left_height"].GetDouble()
        self.g = self.model_part.ProcessInfo[KM.GRAVITY_Z]

    def Check(self):
        super(DryDamBreakBenchmark, self).Check()
        label = "DryDamBreakBenchmark. "
        if self.g <= 0:
            msg = label + "Gravity must be a positive value. Please, check the definition of GRAVITY_Z component in the ProcessInfo."
            raise Exception(msg)
        elif self.hl <= 0:
            msg = label + "Left height must be a positive value. Please, check the Parameters."
            raise Exception(msg)
        elif self.dam <= 0:
            msg = label + "The dam position must be a positive value. Please, check the Parameters."
            raise Exception(msg)

    def Topography(self, coordinates):
        return 0.0

    def Height(self, coordinates, time):
        x = coordinates.X

        xa = self.__xa(time)
        xb = self.__xb(time)

        if x < xa:
            return self.hl
        elif x < xb:
            return 4 / 9 / self.g * (np.sqrt(self.g * self.hl) - 0.5*(x - self.dam) / time)**2
        else:
            return 0.0

    def Velocity(self, coordinates, time):
        x = coordinates.X

        xa = self.__xa(time)
        xb = self.__xb(time)

        if x < xa:
            return [0.0, 0.0, 0.0]
        elif x < xb:
            return [2 / 3 * ((x - self.dam) / time + np.sqrt(self.g * self.hl)), 0.0, 0.0]
        else:
            return [0.0, 0.0, 0.0]

    def __xa(self, t):
        return self.dam - t * np.sqrt(self.g * self.hl)

    def __xb(self, t):
        return self.dam + 2 * t * np.sqrt(self.g * self.hl)
