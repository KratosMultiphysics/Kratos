import KratosMultiphysics as KM

from KratosMultiphysics.ShallowWaterApplication.benchmarks.base_benchmark_process import BaseBenchmarkProcess

# Other imports
from math import sqrt
import scipy.optimize as opt

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DamBreakBenchmark(model, settings["Parameters"])

class DamBreakBenchmark(BaseBenchmarkProcess):
    """Dam break benchark.

    O. Delestre, C. Lucas, P.-A. Ksinant, F. Darboux, C. Laguerre, T.N.T. Vo, F. James, S. Cordier
    SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies
    International Journal for Numerical Methods in Fluids, Wiley, 2013, 72 (3), pp.269-300
    """

    def __init__(self, model, settings ):
        """Constructor of the benchmark.

        The base class validates the settings and sets the model_part, the variables and the benchmark_settings
        """

        super().__init__(model, settings)

        self.dam = self.settings["benchmark_settings"]["dam_position"].GetDouble()
        self.hl = self.settings["benchmark_settings"]["left_height"].GetDouble()
        self.hr = self.settings["benchmark_settings"]["right_height"].GetDouble()
        self.g = self.model_part.ProcessInfo[KM.GRAVITY_Z]

        self.cm = self.__cm()


    def Check(self):
        """This method checks if the input values have physical sense."""

        super().Check()
        label = "DamBreakBenchmark. "
        if self.g <= 0:
            msg = label + "Gravity must be a positive value. Please, check the definition of GRAVITY_Z component in the ProcessInfo."
            raise Exception(msg)
        elif self.hr < 0:
            msg = label + "Right height must be non-negative. Please, check the Parameters."
            raise Exception(msg)
        elif self.hl < 0:
            msg = label + "Left height must be non-negative. Please, check the Parameters."
            raise Exception(msg)


    @classmethod
    def _GetBenchmarkDefaultSettings(cls):
        return KM.Parameters("""
            {
                "dam_position"  : 5.0,
                "left_height"   : 2.0,
                "right_height"  : 1.0
            }
            """
            )


    def _Height(self, coordinates, time):
        x = coordinates.X

        xa = self.__xa(time)
        xb = self.__xb(time)
        xc = self.__xc(time)

        if x < xa:
            return self.hl
        elif x < xb:
            return 4 / 9 / self.g * (sqrt(self.g * self.hl) - 0.5*(x - self.dam) / time)**2
        elif x < xc:
            return self.cm**2 / self.g
        else:
            return self.hr


    def _Velocity(self, coordinates, time):
        x = coordinates.X

        xa = self.__xa(time)
        xb = self.__xb(time)
        xc = self.__xc(time)

        if x < xa:
            return [0.0, 0.0, 0.0]
        elif x < xb:
            return [2 / 3 * ((x - self.dam) / time + sqrt(self.g * self.hl)), 0.0, 0.0]
        elif x < xc:
            return [2 * (sqrt(self.g * self.hl) - self.cm), 0.0, 0.0]
        else:
            return [0.0, 0.0, 0.0]


    def __xa(self, t):
        return self.dam - t * sqrt(self.g * self.hl)


    def __xb(self, t):
        return self.dam + t * (2*sqrt(self.g * self.hl) - 3 * self.cm)


    def __xc(self, t):
        if self.hr > 0:
            return self.dam + t * 2 * self.cm**2 * (sqrt(self.g*self.hl) - self.cm) / (self.cm**2 - self.g * self.hr)
        else:
            return self.__xb(t)


    def __cm(self):
        if self.hr > 0:
            cm0 = sqrt(self.g * 0.5 * (self.hl + self.hr))
            cm = opt.newton(self.__cm_residual, cm0)
            return cm
        else:
            return 0.0


    def __cm_residual(self,cm):
        hl = self.hl
        hr = self.hr
        g  = self.g
        return -8*g*hr*cm**2*(sqrt(g*hl)-cm)**2 + (cm**2-g*hr)**2 * (cm**2+g*hr)
