# Importing the Kratos Library
import KratosMultiphysics as KM

from KratosMultiphysics.ShallowWaterApplication.benchmarks.base_benchmark_process import BaseBenchmarkProcess

# Other imports
import numpy as np
import scipy.optimize as opt

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DamBreakBenchmark(model, settings["Parameters"])

class DamBreakBenchmark(BaseBenchmarkProcess):

    def __init__(self, model, settings ):
        super(DamBreakBenchmark, self).__init__(model, settings)

        benchmark_default_settings = KM.Parameters("""
            {
                "dam_position"  : 5.0,
                "left_height"   : 2.0,
                "right_height"  : 1.0
            }
            """
            )

        self.benchmark_settings.ValidateAndAssignDefaults(benchmark_default_settings)

        self.dam = self.benchmark_settings["dam_position"].GetDouble()
        self.hl = self.benchmark_settings["left_height"].GetDouble()
        self.hr = self.benchmark_settings["right_height"].GetDouble()
        self.g = self.model_part.ProcessInfo[KM.GRAVITY_Z]

        self.cm = self.__cm()

    def Check(self):
        super(DamBreakBenchmark, self).Check()
        label = "DamBreakBenchmark. "
        if self.g <= 0:
            msg = label + "Gravity must be a positive value"
            raise Exception(msg)
        elif self.hr <= 0:
            msg = label + "Right height must be a positive value"
            raise Exception(msg)
        elif self.hl <= 0:
            msg = label + "Left height must be a positive value"
            raise Exception(msg)
        elif self.dam <= 0:
            msg = label + "The dam position must be a positive value"
            raise Exception(msg)

    def Height(self, coordinates, time):
        x = coordinates.X

        xa = self.__xa(time)
        xb = self.__xb(time)
        xc = self.__xc(time)

        if x < xa:
            return self.hl
        elif x < xb:
            return 4 / 9 / self.g * (np.sqrt(self.g * self.hl) - 0.5*(x - self.dam) / time)**2
        elif x < xc:
            return self.cm**2 / self.g
        else:
            return self.hr

    def Velocity(self, coordinates, time):
        x = coordinates.X

        xa = self.__xa(time)
        xb = self.__xb(time)
        xc = self.__xc(time)

        if x < xa:
            return [0.0, 0.0, 0.0]
        elif x < xb:
            return [2 / 3 * ((x - self.dam) / time + np.sqrt(self.g * self.hl)), 0.0, 0.0]
        elif x < xc:
            return [2 * (np.sqrt(self.g * self.hl) - self.cm), 0.0, 0.0]
        else:
            return [0.0, 0.0, 0.0]

    def Momentum(self, coordinates, time):
        return self.Height(coordinates, time) * self.Velocity(coordinates, time)

    def __xa(self, t):
        return self.dam - t * np.sqrt(self.g * self.hl)

    def __xb(self, t):
        return self.dam + t * (2*np.sqrt(self.g * self.hl) - 3 * self.cm)

    def __xc(self, t):
        return self.dam + t * 2 * self.cm**2 * (np.sqrt(self.g*self.hl) - self.cm) / (self.cm**2 - self.g * self.hr)

    def __cm(self):
        cm0 = np.sqrt(self.g * 0.5 * (self.hl + self.hr))
        cm = opt.newton(self.__cm_residual, cm0)
        return cm

    def __cm_residual(self,cm):
        hl = self.hl
        hr = self.hr
        g  = self.g
        return -8*g*hr*cm**2*(np.sqrt(g*hl)-cm)**2 + (cm**2-g*hr)**2 * (cm**2+g*hr)
