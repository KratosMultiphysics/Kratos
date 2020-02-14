# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

from KratosMultiphysics.ShallowWaterApplication.benchmarks.base_benchmark_process import BaseBenchmarkProcess

# Other imports
from numpy import np

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DamBreakBenchmark(model, settings["Parameters"])

class DamBreakBenchmark(BaseBenchmarkProcess):

    def __init__(self, model, settings ):
        super(DamBreakBenchmark, self).__init__()

        default_settings = KM.Parameters("""
            {
                "dam_break_position"  : 5.0,
                "left_height"         : 2.0,
                "right_height"        : 1.0
            }
            """
            )

        self.benchmark_settings.ValidateAndAssignDefaults(default_settings)

        self.dam = self.benchmark_settings["dam_break_position"].GetDouble()
        self.hl = self.benchmark_settings["left_height"].GetDouble()
        self.hr = self.benchmark_settings["right_height"].GetDouble()
        self.g = self.model_part.ProcessInfo[KM.GRAVITY_Z]

        self.cm = __cm()

    def Height(self, coordinates, time):
        x = coordinates.X()

        xa = self.__xa(x, time)
        xb = self.__xb(x, time)
        xc = self.__xc(x, time)

        if x < xa:
            return self.hl
        else if x < xb:
            return 4 / 9 / self.g * (np.sqrt(self.g * self.hl) - 0.5*(x - self.dam) / time)**2
        else if x < xc:
            return self.cm**2 / self.g
        else:
            return self.hr

    def __xa(self, t):
        return self.dam - t * np.sqrt(self.g * self.hl)

    def __xb(self, t):
        return self.dam + t * (2*np.sqrt(self.g * self.hl) - 3 * self.cm)

    def __xc(self, t):
        return self.dam + t * 2 * self.cm**2 * (np.sqrt(self.g*self.hl) - self.cm) / (self.cm**2 - self.g * self.hr)

    def __cm(self):
        raise Exception("TODO: implement this")
        p = np.empty(7)
        p[0] = 1
        p[1] = 0
        p[2] = -9 * self.g * self.hr
        p[3] = 16 * self.g * self.hr * np.sqrt(self.g * self.hl)
        p[4] = -8 * self.g**2 * self.hl * self.hr - self.g**2 * self.hr**2
        p[5] = 0
        p[6] = self.g**3 * self.hr**3
        roots = np.roots(p)
        return roots[0]
